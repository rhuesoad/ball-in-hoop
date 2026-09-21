
import time 
import odrive 
from odrive.enums import(
    AxisState,
    ProcedureResult,
    MotorType, 
    EncoderId,
)
from odrive.utils import dump_errors

# 0) Connection
# Connection avec USB via USB isolator pour trouver l'odrive 

print("Connection to ODrive")
odrv = odrive.find_any()
print("Connected!")
odrv.clear_errors()

# 1) Power-Source limits 
odrv.config.dc_bus_overvoltage_trip_level = 30.0 # 51.0
odrv.config.dc_bus_undervoltage_trip_level = 10.0
odrv.config.dc_max_positive_current = 8.0 # 20.0
odrv.config.dc_max_negative_current = -1.0

# 2) Motor parameters (D5065, 270KV, 7 pole pairs)
odrv.axis0.config.motor.motor_type = MotorType.HIGH_CURRENT
odrv.axis0.config.motor.pole_pairs = 7  # a def pour conversion angle mec / elec
odrv.axis0.config.motor.torque_constant = 8.27/270.0 # --> cf. sheet moteur
# Sets limits of current --> Quite weak to start, but strong enough to overcome friction
#odrv.axis0.config.motor.current_soft_max = 10 
#odrv.axis0.config.motor.current_hard_max = 18

# Idee: ODrive injecte un courant DC pour avoir R = V / I, puis une impulsion pour mesurer L. 
odrv.axis0.config.motor.calibration_current = 6.0 # 10.0 
odrv.axis0.config.motor.resistance_calib_max_voltage = 4.0 

odrv.axis0.config.calibration_lockin.current = 4.5 # 15.0 # Lock-in current quand on mesure l'offset encodeur, max 75% du DC current

# 3) Encoder parameters (AMT102-V, 8192CPR quad incremental)
odrv.inc_encoder0.config.cpr = 8192
odrv.inc_encoder0.config.enabled = True
# Set the communication as commutation encoder (for driver itself) and as load encoder (for feedback)
odrv.axis0.config.load_encoder = EncoderId.INC_ENCODER0
odrv.axis0.config.commutation_encoder = EncoderId.INC_ENCODER0

# RUN A STATE & CHECK RESULT + BLOCK UNTIL IT FINISHES
def run_state(axis, state, name, timeout=30.0):
    "Request a calibration AxisState and block until it finishes"

    # 1) Envoie l'ordre de se mettre ds un state a l'ODrive
    print(f"\n--> {name}")
    axis.requested_state = state            
    t0 = time.time()

    # 2) Attend la fin du process. Si timeout, on a attendu trop longtemps donc on arrete pouru eviter boucle infinie 
    while axis.current_state != AxisState.IDLE:
        if time.time() - t0 > timeout:
            print(f" TIMEOUT waiting for {name}")
            break
        time.sleep(0.2)

    # 3) Verifie le resultat du process 
    if axis.procedure_result == ProcedureResult.SUCCESS:
        print(f"    {name}: SUCCESS")
        return True
    print (f"   {name}: FAILED -> {axis.procedure_result}")
    dump_errors(odrv)
    return False 

# --------------
# 4) MOTOR CALIBRATION
# --------------
# Measure R and L to get iq --> tau = Kt iq
motor_ok = run_state(odrv.axis0, AxisState.MOTOR_CALIBRATION, "MOTOR_CALIBRATION")

if motor_ok:
    R = odrv.axis0.config.motor.phase_resistance
    L = odrv.axis0.config.motor.phase_inductance
    print(f"    Phase Resistance: R = {R*1e3:.3f} mOhm")
    print(f"    Phase Inductance: L = {L*1e6:.3f} uH")
 
# --------------
# 5) ENCODER CALIBRATION
# --------------
if motor_ok:
    input("\n Press Enter when the motor can spin freely ...")
    run_state(odrv.axis0, AxisState.ENCODER_OFFSET_CALIBRATION, "ENCODER_OFFSET_CALIBRATION")

else: 
    print("Encoder calibration skipped because motor calibration failed")
    encoder_ok = False 

# 6) Persist configuration
print("\nSaving Configuration (the ODrive will reboot)")
try: 
    odrv.save_configuration()
except Exception: 
    # save drops the USB link on reboot
    pass
print("Done. Reconnect to verify, then test closed-loop control")
