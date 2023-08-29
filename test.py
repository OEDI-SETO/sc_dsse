
import logging
import helics as h
import json
import numpy as np
from pydantic import BaseModel
from enum import Enum
from typing import List, Optional, Union
from scipy.optimize import least_squares
import scipy.sparse
# import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat 
from scipy.optimize import minimize
import random, sys 
import pickle
from oct2py import octave 

from datetime import datetime
from oedisi.types.data_types import (
    AdmittanceSparse,
    MeasurementArray,
    AdmittanceMatrix,
    Topology,
    Complex,
    VoltagesMagnitude,
    VoltagesAngle,
    VoltagesReal,
    VoltagesImaginary,
    PowersReal,
    PowersImaginary,
)


logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


def power_flow_equations(x, Y_bus):
    num_buses = Y_bus.shape[0]
    V_magnitude = x[:num_buses]
    V_angle = x[num_buses:]

    V_complex = V_magnitude * np.exp(1j * V_angle)
    S = V_complex * np.conj(Y_bus @ V_complex)

    P = S.real
    Q = S.imag


    return P, Q

def cal_h(knownP, knownQ, knownV, Y, deltaK, VabsK, num_node):
    h1 = (VabsK[knownV]).reshape(-1,1)
    Vp = VabsK * np.exp(1j * deltaK)
    S = Vp * (Y.conjugate() @ Vp.conjugate())
    P, Q = S.real, S.imag
    h2, h3 = P[knownP].reshape(-1,1), Q[knownQ].reshape(-1,1)
    h = np.concatenate((h1, h2, h3), axis=0)
    return h.reshape(-1)

def cal_H(X0, z, num_node, knownP, knownQ, knownV, Y):
    deltaK, VabsK = X0[:num_node], X0[num_node:]
    num_knownV = len(knownV)
    #Calculate original H1
    H11, H12 = np.zeros((num_knownV, num_node)), np.zeros(num_knownV * num_node)
    H12[np.arange(num_knownV)*num_node + knownV] = 1
    H1 = np.concatenate((H11, H12.reshape(num_knownV, num_node)), axis=1)
    Vp = VabsK * np.exp(1j * deltaK)
##### S = np.diag(Vp) @ Y.conjugate() @ Vp.conjugate()
######  Take gradient with respect to V
    H_pow2 = (
            Vp.reshape(-1, 1) * Y.conjugate() * np.exp(-1j * deltaK).reshape(1, -1) + 
        np.exp(1j * deltaK) * np.diag(Y.conjugate() @ Vp.conjugate())
        )
    # Take gradient with respect to delta
    H_pow1 = (
            1j * Vp.reshape(-1, 1) * (np.diag(Y.conjugate() @ Vp.conjugate()) -
                Y.conjugate() * Vp.conjugate().reshape(1, -1))
            )
        
    H2 = np.concatenate((H_pow1.real, H_pow2.real), axis=1)[knownP, :]
    H3 = np.concatenate((H_pow1.imag, H_pow2.imag), axis=1)[knownQ, :]
    H = np.concatenate((H1, H2, H3), axis=0)   
    return -H

def residual(X0, z, num_node, knownP, knownQ, knownV, Y):
    delta, Vabs = X0[:num_node], X0[num_node:]
    h = cal_h(knownP, knownQ, knownV, Y, delta, Vabs, num_node)
    return z-h
def residual2(X0, z, num_node, knownP, knownQ, knownV, Y):
    y = residual(X0, z, num_node, knownP, knownQ, knownV, Y)
    y2 = sum(y**2)
    return y2



def jac2(X0, z, num_node, knownP, knownQ, knownV, Y):
    residual_cal = residual(X0, z, num_node, knownP, knownQ, knownV, Y)
    jac_cal = cal_H(X0, z, num_node, knownP, knownQ, knownV, Y)
    jac2 = np.transpose(residual_cal)@jac_cal
    jac3 = 2*jac2
    return jac3


def matrix_to_numpy(admittance: List[List[Complex]]):
    "Convert list of list of our Complex type into a numpy matrix"
    return np.array([[x[0] + 1j * x[1] for x in row] for row in admittance])


class UnitSystem(str, Enum):
    SI = "SI"
    PER_UNIT = "PER_UNIT"


class Power:
    def __init__(self, ids, values):
        self.ids = ids
        self.values = values


class AlgorithmParameters(BaseModel):
    tol: float = 5e-7
    units: UnitSystem = UnitSystem.PER_UNIT
    base_power: Optional[float] = 100.0

    class Config:
        use_enum_values = True


def get_indices(topology, measurement):
    "Get list of indices in the topology for each index of the input measurement"
    inv_map = {v: i for i, v in enumerate(topology.base_voltage_magnitudes.ids)}
    return [inv_map[v] for v in measurement.ids]

# def get_indices_aux(topology, measurement):
#     "Get list of indices in the topology for each index of the input measurement"
#     inv_map = {v: i for i, v in enumerate(topology.base_voltage_magnitudes.ids)}
#     return [inv_map[v] for v in measurement.ids]


def state_estimator_dsse(P_MEAS, Q_MEAS, VMAG_MEAS, KNOWNP, KNOWNQ, KNOWNV, YBUS):
    '''
    ALL VALUES ARE PER UNIT
         
    P_MEAS: Active power (per unit) 
    Q_MEAS: Reactive power (per unit)
    VMAG_MEAS: Voltage Magnitude (per unit)
    KNOWNP, KNOWNQ, KNOWNV: Indices of known measurements
    YBUS: Ybus (per unit)
    '''
    # print ("types P,Q,V,KP,Y: ", type(P_MEAS), type(Q_MEAS),type(VMAG_MEAS),type(KNOWNP), type(YBUS)) 



    Z = np.concatenate((VMAG_MEAS, P_MEAS, Q_MEAS), axis=0)
    num_node_orig = len(YBUS)
    # initialize 
    initial_ang=0
    initial_V=1
    delta = np.full(num_node_orig, initial_ang)
    Vabs = np.full(num_node_orig , initial_V)
    X0 = np.concatenate((delta, Vabs))

    low_limit = np.concatenate((np.ones(num_node_orig)* (- np.pi - np.pi/6),
                            np.ones(num_node_orig)*0.95))
    up_limit = np.concatenate((np.ones(num_node_orig)* (np.pi + np.pi/6),
                                np.ones(num_node_orig)*1.05))

    # Modify bounds for node 1
    low_limit[0] = 0  # Angle should be 0
    up_limit[0] = 0.000001
    low_limit[num_node_orig] = 1  # Magnitude should be 1
    up_limit[num_node_orig] = 1.000001

    bounds_1 = tuple(zip(low_limit,up_limit))

    tol= 1e-10

    res_1 = minimize(
        residual2,
        X0,
        jac=jac2,
        bounds=bounds_1,
        # method = 'SLSQP',
        # method='L-BFGS-B',
        # method = 'Newton-CG',
        method = 'BFGS',
        # CG was also working but slow
        options={'disp': True,'gtol': tol, 'ftol': tol, 'maxiter':1000},
        args=(Z, num_node_orig, KNOWNP, KNOWNQ, KNOWNV, YBUS),
    )
    # Few comments: 
    # 1. Both BFGS and CG are giving good results. BFGS is faster. 
    # 2. BFGS doesn't support bounds but still works. note that the choice of  Sbase chnages resutls. 
    # 3. SLSQP supports bounds and its options={'disp': True,'gtol': tol, 'ftol': tol, 'maxiter':500}, works
    # 4. SLSQP is slow 
    # 5. L-BFGS-B can also be made to give better results. Needs tuning by selecting good params. (read documentation)



    result1 = res_1.x
    vmagestDecen1, vangestDecen1 = result1[num_node_orig:], result1[:num_node_orig]
    return vmagestDecen1, vangestDecen1


def get_y(admittance: Union[AdmittanceMatrix, AdmittanceSparse], ids: List[str]):
    if type(admittance) == AdmittanceMatrix:
        assert ids == admittance.ids
        return matrix_to_numpy(admittance.admittance_matrix)
    elif type(admittance) == AdmittanceSparse:
        node_map = {name: i for (i, name) in enumerate(ids)}
        return scipy.sparse.coo_matrix(
            (
                [v[0] + 1j * v[1] for v in admittance.admittance_list],
                (
                    [node_map[r] for r in admittance.from_equipment],
                    [node_map[c] for c in admittance.to_equipment],
                ),
            )
        ).toarray()


class StateEstimatorFederate:
    "State estimator federate. Wraps state_estimation with pubs and subs"

    def __init__(
        self, federate_name, algorithm_parameters: AlgorithmParameters, input_mapping
    ):
        "Initializes federate with name and remaps input into subscriptions"
        deltat = 0.1

        self.algorithm_parameters = algorithm_parameters

        # Create Federate Info object that describes the federate properties #
        fedinfo = h.helicsCreateFederateInfo()

        fedinfo.core_name = federate_name
        fedinfo.core_type = h.HELICS_CORE_TYPE_ZMQ
        fedinfo.core_init = "--federates=1"
        h.helicsFederateInfoSetTimeProperty(
            fedinfo, h.helics_property_time_delta, deltat
        )

        self.vfed = h.helicsCreateValueFederate(federate_name, fedinfo)
        logger.info("Value federate created")

        # Register the publication #
        self.sub_voltages_magnitude = self.vfed.register_subscription(
            input_mapping["voltages_magnitude"], "V"
        )
        self.sub_power_P = self.vfed.register_subscription(
            input_mapping["powers_real"], "W"
        )
        self.sub_power_Q = self.vfed.register_subscription(
            input_mapping["powers_imaginary"], "W"
        )
        self.sub_topology = self.vfed.register_subscription(
            input_mapping["topology"], ""
        )
        self.pub_voltage_mag = self.vfed.register_publication(
            "voltage_mag", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_voltage_angle = self.vfed.register_publication(
            "voltage_angle", h.HELICS_DATA_TYPE_STRING, ""
        )

    def run(self):
        "Enter execution and exchange data"
        # Enter execution mode #
        self.vfed.enter_executing_mode()
        logger.info("Entering execution mode")

        granted_time = h.helicsFederateRequestTime(self.vfed, h.HELICS_TIME_MAXTIME)

        self.initial_ang = None
        self.initial_V = None
        topology = Topology.parse_obj(self.sub_topology.json)
        ids = topology.base_voltage_magnitudes.ids
        logger.info("Topology has been read")
        slack_index = None
        if not isinstance(topology.admittance, AdmittanceMatrix) and not isinstance(
            topology.admittance, AdmittanceSparse
        ):
            raise "Weighted Least Squares algorithm expects AdmittanceMatrix/Sparse as input"

        for i in range(len(ids)):
            if ids[i] == topology.slack_bus[0]:
                slack_index = i

        while granted_time < h.HELICS_TIME_MAXTIME:

            if not self.sub_voltages_magnitude.is_updated():
                granted_time = h.helicsFederateRequestTime(
                    self.vfed, h.HELICS_TIME_MAXTIME
                )
                continue

            logger.info("start time: " + str(datetime.now()))

            voltages = VoltagesMagnitude.parse_obj(self.sub_voltages_magnitude.json)
            power_P = PowersReal.parse_obj(self.sub_power_P.json)
            power_Q = PowersImaginary.parse_obj(self.sub_power_Q.json)
            
            
            # print("voltages:", voltages)
            # print("Power P:", power_P)
            # print("Power Q:", power_Q)
            
            
            # ajay  
            # aux_nodes = [ '149.1', '149.2', '149.3', '1.2', '1.3', '3.3', '8.1', '8.2', '8.3', 
                # '13.1', '13.2', '13.3', '14.1', '18.1', '18.2', '18.3', '15.3', '21.1', 
                # '21.2', '21.3', '23.1', '23.2', '23.3', '25.1', '25.2', '25.3', '26.1', '26.3', 
                # '28.2', '28.3', '27.1', '27.3', '30.1', '30.2', '250.1', '250.2', '250.3', '36.1', '36.2', 
                # '35.3', '40.1', '40.2', '40.3', '42.2', '42.3', '44.1', '44.2', '44.3', '50.1', '50.2', '51.2', '51.3', 
                # '151.1', '151.2', '151.3', '52.2', '52.3', '53.2', '53.3', '54.1', '54.2', '54.3', '57.1', '57.2', '57.3', 
                # '56.1', '56.3', '60.2', '60.3', '61.1', '61.2', '61.3', '62.1', '62.2', '64.1', '64.3', '66.1', '66.2', '67.1', 
                # '67.2', '67.3', '72.1', '72.2', '72.3', '97.1', '97.2', '97.3', '77.1', '77.3', '86.1', '86.3', '78.1', '78.2', 
                # '78.3', '79.2', '79.3', '81.1', '81.2', '81.3', '82.2', '82.3', '83.1', '83.2', '89.1', '89.2', '89.3', '91.1', 
                # '91.2', '91.3', '93.1', '93.2', '93.3', '95.1', '95.3', '98.2', '98.3', '99.1', '99.3', '100.1', '100.2', '450.1', 
                # '450.2', '450.3', '197.1', '197.2', '197.3', '101.1', '101.2', '101.3', '105.1', '105.2', '105.3', '108.1', '108.2', 
                # '108.3', '300.1', '300.2', '300.3', '110.1', '135.1', '135.2', '135.3', '152.1', '152.2', '152.3', '160.1', '160.2', 
                # '160.3', '610.1', '610.2', '610.3']

            # power_P_aux = Power(ids=aux_nodes, values=[0]*len(aux_nodes))
            # power_Q_aux = power_P_aux

            # power_P.ids.extend(power_P_aux.ids)
            # power_P.values.extend(power_P_aux.values)

            # power_Q.ids.extend(power_Q_aux.ids)
            # power_Q.values.extend(power_Q_aux.values)
            
            
            
            pseudo_node_names =  ['65.2', '65.3', '76.2', '76.3']
            pseudo_values_P = [0.72*35/100, 0.74*70/100, 0.73*70/100, 0.76*70/100]
            pseudo_values_Q = [0.72*25/100, 0.74*50/100, 0.73*50/100, 0.76*50/100]

            pseudo_node_P = Power(ids=pseudo_node_names, values=pseudo_values_P)
            pseudo_node_Q = Power(ids=pseudo_node_names, values=pseudo_values_Q)



            remaining_node_names =  ['150R.1', '150R.2', '150R.3', '149.1', '149.2', '149.3', '1.2', '1.3', '3.3', '8.1', '8.2', '8.3', 
                                     '13.1', '13.2', '13.3', '9R.1', '14.1', '18.1', '18.2', '18.3', '15.3', '21.1', '21.2', '21.3', '23.1', 
                                     '23.2', '23.3', '25.1', '25.2', '25.3', '25R.1', '25R.3', '26.1', '26.3', '28.2', '28.3', '27.1', '27.3', 
                                     '30.1', '30.2', '250.1', '250.2', '250.3', '36.1', '36.2', '35.3', '40.1', '40.2', '40.3', '42.2', '42.3', 
                                     '44.1', '44.2', '44.3', '50.1', '50.2', '51.2', '51.3', '151.1', '151.2', '151.3', '52.2', '52.3', '53.2', '53.3', 
                                     '54.1', '54.2', '54.3', '57.1', '57.2', '57.3', '56.1', '56.3', '60.2', '60.3', '61.1', '61.2', '61.3', '62.1', '62.2', 
                                     '64.1', '64.3', '66.1', '66.2', '67.1', '67.2', '67.3', '72.1', '72.2', '72.3', '97.1', 
                                    '97.2', '97.3', '77.1', '77.3', '86.1', '86.3', '78.1', '78.2', '78.3', '79.2', '79.3', '81.1', '81.2', '81.3', '82.2', 
                                    '82.3', '83.1', '83.2', '89.1', '89.2', '89.3', '91.1', '91.2', '91.3', '93.1', '93.2', '93.3', '95.1', '95.3', '98.2', 
                                    '98.3', '99.1', '99.3', '100.1', '100.2', '450.1', '450.2', '450.3', '197.1', '197.2', '197.3', '101.1', '101.2', '101.3', 
                                    '105.1', '105.2', '105.3', '108.1', '108.2', '108.3', '300.1', '300.2', '300.3', '110.1', '135.1', '135.2', '135.3', '152.1', 
                                    '152.2', '152.3', '160R.1', '160R.2', '160R.3', '160.1', '160.2', '160.3', '61S.1', '61S.2', '61S.3', '300_OPEN.1', '300_OPEN.2', 
                                    '300_OPEN.3', '94_OPEN.1', '610.1', '610.2', '610.3']

            
            remaining_node_names = []

            remaining_node_P = Power(ids=remaining_node_names, values=[0.0]*len(remaining_node_names))
            remaining_node_Q = Power(ids=remaining_node_names, values=[0.0]*len(remaining_node_names))
            
            power_P.ids.extend(pseudo_node_P.ids)
            power_P.values.extend(pseudo_node_P.values)

            power_Q.ids.extend(pseudo_node_Q.ids)
            power_Q.values.extend(pseudo_node_Q.values)
            
            
            power_P.ids.extend(remaining_node_P.ids)
            power_P.values.extend(remaining_node_P.values)

            power_Q.ids.extend(remaining_node_Q.ids)
            power_Q.values.extend(remaining_node_Q.values)
            
            
            knownP = get_indices(topology, power_P)
            knownQ = get_indices(topology, power_Q)
            knownV = get_indices(topology, voltages)


            if self.initial_V is None:
                # Flat start or using average measurements
                if len(knownP) + len(knownV) + len(knownQ) > len(ids) * 2:
                    self.initial_V = 1.0
                else:
                    self.initial_V = np.mean(
                        np.array(voltages.values)
                        / np.array(topology.base_voltage_magnitudes.values)[knownV]
                    )
            if self.initial_ang is None:
                self.initial_ang = np.array(topology.base_voltage_angles.values)


            base_voltages = np.array(topology.base_voltage_magnitudes.values)
            ids = topology.base_voltage_magnitudes.ids

            base_power = 100
            
            if parameters.base_power != None:
                base_power = parameters.base_power
# ajay : python based DSSE


            
            Y = get_y(topology.admittance, ids)
            # Hand-crafted unit conversion (check it, it works)
            Y = (
                base_voltages.reshape(1, -1)
                * Y
                * base_voltages.reshape(-1, 1)
                / (base_power * 1000)
            )

            
            # ajay another thought
            # size_index = len(Y)

            # total_index = list(range(size_index))

            # # convert the list to set
            # total_index_set = set(total_index)
            # knownP_set = set(knownP)

            # not_in_knownP = total_index_set - knownP_set
            # power_P_values = power_P.values 

            # for index in not_in_knownP:
            #     knownP.append(index)
            #     power_P_values.append(0)
            
            # # DO THE SAME FOR Q   
                                 
            # # total_index_set = set(total_index)
            # knownQ_set = set(knownQ)

            # not_in_knownQ = total_index_set - knownQ_set
            # power_Q_values = power_Q.values 

            # for index in not_in_knownQ:
            #     knownQ.append(index)
            #     power_Q_values.append(0)



            z = np.concatenate(
                (
                    voltages.values / base_voltages[knownV],
                    -np.array(power_P.values) / base_power,
                    -np.array(power_Q.values) / base_power,
                ),
                axis=0,
            )    
            # as per the construction of z in the previous line 
            power_P_values = -np.array(power_P.values) / base_power
            power_Q_values = -np.array(power_Q.values) / base_power

            voltage_values = voltages.values / base_voltages[knownV]

            # print('Topology IDS: ', topology.base_voltage_magnitudes.ids)
            # print('Entire Topology: ', topology)
            # print('IDS from Ybus: ', topology.base_voltage_magnitudes.ids)


            vars_to_save = {"power_P": power_P.values, "power_Q": power_Q.values, "voltage": voltages.values, 
            "Y": Y,"ids_PQ": power_P.ids, "base_power": base_power, "base_voltages": base_voltages, 
            "ids_voltage": voltages.ids, "knownP": knownP, "knownQ": knownQ, "knownV": knownV, "NodeIDs":ids}
            
            savemat('saved_vars.mat', vars_to_save)
            
            # with open('saved_vars.pkl', 'wb') as f:
            #     pickle.dump(vars_to_save, f)
            
            # voltage_values = []
            # knownV = []
            # AJAY : note that values must be in PU
            
            # voltage_magnitudes, voltage_angles =\
            #     state_estimator_dsse(power_P_values, power_Q_values, voltage_values, knownP, knownQ, knownV, Y)
            
            #     # state_estimator_dsse(power_P.values, power_Q.values, voltages.values, knownP, knownQ, knownV, Y)
            #     #state_estimator_dsse(P_MEAS, Q_MEAS, VMAG_MEAS, KNOWNP, KNOWNQ, KNOWNV, YBUS)

            # print("voltage magnitudes: ", voltage_magnitudes)    

            # JENNY CODE
            # extract the variables
            # NodeIDs = vars_to_save['NodeIDs']
            # ids_PQ = vars_to_save['ids_PQ']
            # ids_voltage = vars_to_save['ids_voltage']
            # power_P = vars_to_save['power_P']
            # power_Q = vars_to_save['power_Q']
            # voltage = vars_to_save['voltage']
            # Y = vars_to_save['Y']
            # base_power = vars_to_save['base_power']

            # output = octave.jennyPreprocess(NodeIDs, ids_PQ, ids_voltage, power_P, power_Q, voltage, Y, base_power)
            output = octave.jennyPreprocess2()

            # # # CHECK THIS     
            voltage_magnitudes = output['VN']
            voltage_angles = output['AN']


            # print("voltage mag: ", output['VN'])
            # print("voltage mag: ", type(output['VN']))
            # print("type of voltage angles: ", type(voltage_angles))


            #voltage_magnitudes, voltage_angles = state_estimator(
            #    self.algorithm_parameters,
            #    topology,
            #    power_P,
            #    power_Q,
            #    voltages,
            #    initial_V=self.initial_V,
            #    initial_ang=self.initial_ang,
            #    slack_index=slack_index,
            #)

            # self.initial_V = voltage_magnitudes
            # self.initial_ang = voltage_angles
            self.pub_voltage_mag.publish(
                VoltagesMagnitude(
                    values=list(voltage_magnitudes), ids=ids, time=voltages.time
                ).json()
            )
            self.pub_voltage_angle.publish(
                VoltagesAngle(
                    values=list(voltage_angles), ids=ids, time=voltages.time
                ).json()
            )
            logger.info("end time: " + str(datetime.now()))

        self.destroy()

    def destroy(self):
        "Finalize and destroy the federates"
        h.helicsFederateDisconnect(self.vfed)
        logger.info("Federate disconnected")

        h.helicsFederateFree(self.vfed)
        h.helicsCloseLibrary()


if __name__ == "__main__":
    with open("static_inputs.json") as f:
        config = json.load(f)
        federate_name = config["name"]
        if "algorithm_parameters" in config:
            parameters = AlgorithmParameters.parse_obj(config["algorithm_parameters"])
        else:
            parameters = AlgorithmParameters.parse_obj({})

    with open("input_mapping.json") as f:
        input_mapping = json.load(f)

    sfed = StateEstimatorFederate(federate_name, parameters, input_mapping)
    sfed.run()

