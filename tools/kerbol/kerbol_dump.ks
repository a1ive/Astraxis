// Dumps the stock Kerbol system (orbits, physical parameters, rotation,
// atmospheres) for the Astraxis Kerbol scene.
//
// Copy to <KSP>/Ships/Script/, launch any vessel with a kOS processor, open its
// terminal and run:
//     SWITCH TO 0. RUNPATH("kerbol_dump").
// Output: <KSP>/Ships/Script/kerbol_bodies.txt
//
// Format: "key = value" lines, one [[body]] block per body; atm_sample repeats
// (altitude m, pressure atm, temperature K). Units are KSP's: m, s, kg, deg.
// Vectors are in kOS's raw (ship-centered, left-handed, possibly rotating)
// frame; each body block records the solar prime vector and Kerbol's angular
// velocity sampled in the same physics tick as its position, so the frame can
// be resolved offline. sample_ut_start/end bracket that tick (they should match).

@LAZYGLOBAL OFF.
SET CONFIG:IPU TO 2000.

GLOBAL dump_path IS "0:/kerbol_bodies.txt".
IF EXISTS(dump_path) { DELETEPATH(dump_path). }

FUNCTION emit { PARAMETER s. LOG s TO dump_path. }
FUNCTION quoted { PARAMETER s. RETURN CHAR(34) + s + CHAR(34). }
// Names must not shadow kOS built-ins (v, r, q, ...): CLOBBERBUILTINS.
FUNCTION vec_str { PARAMETER vec. RETURN "[" + vec:X + ", " + vec:Y + ", " + vec:Z + "]". }

emit("# kerbol_dump.ks output (Astraxis tools/kerbol)").
emit("kos_version = " + quoted(CORE:VERSION)).
emit("ut = " + TIME:SECONDS).

LOCAL bodies IS LIST().
LIST BODIES IN bodies.
FOR b IN bodies {
    // Sample everything frame-dependent at the start of one physics tick.
    WAIT 0.
    LOCAL t0 IS TIME:SECONDS.
    LOCAL prime IS SOLARPRIMEVECTOR.
    LOCAL sun_spin IS SUN:ANGULARVEL.
    LOCAL spin IS b:ANGULARVEL.
    LOCAL rot IS b:ROTATIONANGLE.
    LOCAL pos IS V(0, 0, 0).
    LOCAL vel IS V(0, 0, 0).
    LOCAL parent_vel IS V(0, 0, 0).
    IF b:HASBODY {
        SET pos TO b:POSITION - b:BODY:POSITION.
        SET vel TO b:VELOCITY:ORBIT.
        SET parent_vel TO b:BODY:VELOCITY:ORBIT.
    }
    LOCAL t1 IS TIME:SECONDS.

    emit("").
    emit("[[body]]").
    emit("name = " + quoted(b:NAME)).
    emit("mass_kg = " + b:MASS).
    emit("mu_m3_s2 = " + b:MU).
    emit("radius_m = " + b:RADIUS).
    emit("rotation_period_s = " + b:ROTATIONPERIOD).
    emit("has_ocean = " + b:HASOCEAN).
    emit("has_solid_surface = " + b:HASSOLIDSURFACE).

    emit("sample_ut_start = " + t0).
    emit("sample_ut_end = " + t1).
    emit("rotation_angle_deg = " + rot).
    emit("solar_prime_vector_raw = " + vec_str(prime)).
    emit("sun_angular_velocity_raw = " + vec_str(sun_spin)).
    emit("angular_velocity_raw = " + vec_str(spin)).

    IF b:HASBODY {
        LOCAL o IS b:ORBIT.
        emit("parent = " + quoted(b:BODY:NAME)).
        emit("soi_radius_m = " + b:SOIRADIUS).
        emit("semi_major_axis_m = " + o:SEMIMAJORAXIS).
        emit("eccentricity = " + o:ECCENTRICITY).
        emit("inclination_deg = " + o:INCLINATION).
        emit("lan_deg = " + o:LONGITUDEOFASCENDINGNODE).
        emit("arg_periapsis_deg = " + o:ARGUMENTOFPERIAPSIS).
        emit("mean_anomaly_at_epoch_deg = " + o:MEANANOMALYATEPOCH).
        emit("epoch_s = " + o:EPOCH).
        emit("period_s = " + o:PERIOD).
        emit("position_rel_parent_raw_m = " + vec_str(pos)).
        emit("velocity_orbit_raw_m_s = " + vec_str(vel)).
        emit("parent_velocity_orbit_raw_m_s = " + vec_str(parent_vel)).
    }

    IF b:ATM:EXISTS {
        LOCAL atm IS b:ATM.
        emit("atm_height_m = " + atm:HEIGHT).
        emit("atm_sea_level_pressure_atm = " + atm:SEALEVELPRESSURE).
        emit("atm_oxygen = " + atm:OXYGEN).
        emit("atm_molar_mass_kg_mol = " + atm:MOLARMASS).
        emit("atm_adiabatic_index = " + atm:ADIABATICINDEX).
        LOCAL n IS 200.
        FROM { LOCAL i IS 0. } UNTIL i > n STEP { SET i TO i + 1. } DO {
            LOCAL h IS atm:HEIGHT * i / n.
            emit("atm_sample = [" + h + ", " + atm:ALTITUDEPRESSURE(h) + ", " + atm:ALTITUDETEMPERATURE(h) + "]").
        }
    }
    PRINT "dumped " + b:NAME.
}
PRINT "wrote " + dump_path.
