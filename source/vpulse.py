###########################################################
#
#
#               Imports
#
#
###########################################################
print 'V_masetr pulse is connected'
from lib.math import fit
import math
import ATS9360.DataTreatment as dt
import numpy as np
import qt
import datetime
import time
import sys
import matplotlib.pyplot as plt

###########################################################
#
#
#               Devices
#
#
###########################################################

Tabor           = qt.instruments.get('Tabor')
smb_cavity     = qt.instruments.get('smb_1')
smb_atom        = qt.instruments.get('smb_2')
ats9360        = qt.instruments.get('ats9360')
SSB_cavity      = qt.instruments.get('SSB_cavity')
SSB_atom      = qt.instruments.get('SSB_atom')
SSB_atom2     = qt.instruments.get('SSB_atom2')
SSB_3          =qt.instruments.get('SSB_3') #not sure !V 180905

Rudat = qt.instruments.get('RUDAT_ph')
Cur_sour = qt.instruments.get('hp3245')

Pulsing_instrument = qt.instruments.get('Pulsing_instrument')

SSB_cavity.set_freq_start(4)
SSB_cavity.set_freq_stop(8)
SSB_cavity.set_conversion_loss(6.)
SSB_cavity.set_LO_power(15)
SSB_cavity.set_band_type(-1)
SSB_cavity.set_IF_frequency(0.0)

SSB_atom.set_freq_start(4)
SSB_atom.set_freq_stop(8.)
SSB_atom.set_conversion_loss(6.)
SSB_atom.set_LO_power(15)
SSB_atom.set_band_type(-1)
SSB_atom.set_IF_frequency(0.05)

SSB_3.set_freq_start(4)
SSB_3.set_freq_stop(8.)
SSB_3.set_conversion_loss(6.)
SSB_3.set_LO_power(15)
SSB_3.set_band_type(-1)
SSB_3.set_IF_frequency(0.05)



###########################################################
#
#
#               Class of parameters
#
#
###########################################################

class SetOfParam:
    #readout tone
    power1 = -0.0
    freq_read = 7.034 #GHz
    t_read = 500e-9 #duration of readout pulse # in sec
    try:
        rudat = Rudat.get_attenuation()
    except:
        rudat = 0
    #qubit 0-1 tone
    power2 = 0
    freq_q = 6.284  #GHz
    tpi = 30e-9 # in sec
    #other
    try:
        cur = Cur_sour.get_current()
    except:
        cur = 0
    nsigma = 0

    def show(self):
        print ' power1   =', self.power1,'\n freq_read=', self.freq_read, 'GHz'
        print  ' t_read   =', 1e9*self.t_read,'ns', '\n Rudat    =', self.rudat,'dB', '\n \n power2   =', self.power2,'\n freq_q   =', self.freq_q,'GHz', '\n tpi      =', 1e9*self.tpi,'ns'
        print ' nsigma   =', self.nsigma
        print '\n current  =', self.cur, 'mA'

    def get_string(self):
        name = ' p1='+str(round(self.power1,2))+';fr='+str(round(self.freq_read,3))+';tr='+str(round(self.t_read,0))
        name = name +';Rd='+str(round(self.rudat,2))+';p2='+str(round(self.power2,2))
        name = name +  ';fq='+str(round(self.freq_q,3))+';pi='+str(round(self.tpi))+';cur'+str(round(self.cur,4))+'sig='+str(round(self.nsigma,1))
        return name
    #def copy(self):
    #    return self

###########################################################
#
#
#               Handmade math functions
#
#
###########################################################

def closest_even(x):
    '''
    Function to find the closest even number to x
    x should be more than 1.0
    '''
    if abs(x) < 1.0:
        print 'error: closest_even() - abs(x) must be more than 1.0'
    if x < 0:
        sign = -1
    else:
        sign = 1
    x = abs(x)

    # if rounds to even - simple case
    if ( round(x) %2 == 0):
        result = round(x)
        return result*sign
    # round(x) is not even:
    #two even candidats:
    a1 = round(x)+1
    a2 = round(x)-1
    r1 = abs(a1-x)
    r2 = abs(a2-x)
    if (r1 < r2):
        result = a1
    else:
        result = a2
    return result*sign

###########################################################
#
#
#               Measurment functions
#
#
###########################################################
def get_rabi_pi(parameters, power2, Tabor_loading = 1, averaging=1e3, nom_expect = 0, window = 400):
    '''
    Return result in seconds: tpi = 30e-9
    '''
    qt.msleep(1)
    FIT = True
    Ajust_window = False

    try:
        GausSigma = parameters.nsigma
        tr_start = 0e-9
        tr_stop = window*1e-9
        tr_step = tr_stop/200.
        power1 = parameters.power1
        #power2 = parameters.power2
        f_atom = parameters.freq_q
        f_cav = parameters.freq_read
        t_meas = parameters.t_read

    except:
        print 'WARNING: (rabi) one of parameters is not defined!'
        return None

    try:
        delta_t = 200
        acq_time = t_meas*1e9 + delta_t + 300.

        #tr_stop = 0.4e-6/2 #increased 2 times - ok
        #tr_step = 1e-9/2
        #tr_start = 0e-9
        t_wait = 0e-9

        T_vec = np.arange(tr_start, tr_stop, tr_step)
        nb_sequences = len(T_vec)

        t_rise = None
    except:
        print 'WARNING: (rabi) problem of setting'
        return None

    ###########################################################
    try:
        mw = 2
        Pulsing_instrument.set_routing_awg({'secondtone_channel':4})

        if Tabor_loading:
            if COMMENTS:
                print 'Tabor loading'
            Pulsing_instrument.set_trigger_time(100.)
            if GausSigma < 0:
                GausSigma = 0
                print 'Warning_Fidelity: Gaussigma was <0, changed to 0'
            if GausSigma > 5:
                GausSigma = 5
                print 'Warning_Fidelity: Gaussigma was >5, changed to 5'
            if GausSigma > 0:
                    if COMMENTS:
                        print 'gaussian form!'
            #gaussian pulse here
            Pulsing_instrument.write_Rabi_pulsessequence(tr_stop, tr_step, tr_start, t_meas,
                t_wait=t_wait, delta_m1_start=0.1e-6, phi=0, delete='all', t_rise =t_rise, nsigma=GausSigma )

        qt.mstart()
    except:
        print 'WARNING: (rabi) problem of Tabor loading'
        return None

    data_measurement = qt.Data(name='Rabi' + parameters.get_string() +'nsigma='+ str(GausSigma))
    data_measurement.add_coordinate('excitation time [ns]', units = 'ns')
    data_measurement.add_value('S21 ',            units = 'Volt')
    data_measurement.add_value('Phase ',            units = 'rad')
    data_measurement.add_value('Re ',            units = 'Volt')
    data_measurement.add_value('Im ',            units = 'Volt')
    data_measurement.create_file()

    try:
        plot2d_1 = qt.Plot2D(data_measurement,
                          name      = 'S21 rabi',
                          coorddim  = 0,
                          valdim    = 1)

        plot2d_2 = qt.Plot2D(data_measurement,
                            name      = 'Phase rabi',
                            coorddim  = 0,
                            valdim    = 2,
                            maxtraces = 2)
        #
        plot2d_3 = qt.Plot2D(data_measurement,
                          name      = 'Re rabi',
                          coorddim  = 0,
                          valdim    = 3,
                          maxtraces = 2)

        plot2d_4 = qt.Plot2D(data_measurement,
                            name      = 'Im rabi',
                            coorddim  = 0,
                            valdim    = 4,
                            maxtraces = 2)
        if FIT:
            data_fit = qt.Data(name='Rabi_OSC_fit')
            data_fit.add_value('parameters ',            units = 'rad, rad, GHz*2pi, rad, ns')
            data_fit.add_value('errors ',            units = 'rad, rad, GHz*2pi, rad, ns')
            data_fit.create_file()
            # TT_vec = np.append(T_vec, T_vec)
        # board_flag = None
    except:
        print 'WARNING: (rabi) problem with plot the graph'

    try:
        Pulsing_instrument.prep_rabi(f_cav, f_atom, averaging, nb_sequences,
            power1, power2, acq_time, t_meas*1e9, delta_t, mw =mw)
        if COMMENTS:
            print 'Rabi_sequence was prepared'
        qt.msleep(2.1)
        Tabor.set_trigger_source('TIM')
        print '2'#!V
        cycle_counter = 0
        while Pulsing_instrument.get_acquisition_completed() !=100.:
            #print  Pulsing_instrument.get_acquisition_completed(), '%' #!V
            result = Pulsing_instrument.measurement()
            ((real_a, rea0), (imag_a, ima0))= result
            real_a -= rea0
            imag_a -= ima0 #deduction of background V?
            amplitude = np.sqrt(real_a**2+imag_a**2)
            complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)

            qt.msleep(0.1)

            plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
            plot2d_2.replace_inline_data(T_vec*1e9, phase)
            # plot2d_3.replace_inline_data(T_vec*1e9, real_a)
            # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)
            #old if FIT
            if False:
                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, amplitude)
                # guess parameters##########################################################
                background = (amplitude.max() + amplitude.min() )/2.
                osc_amp = (amplitude.max() - amplitude.min() )/2.
                nb_expected_oscillation = nom_expect
                pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
                # print pulsation/2/np.pi

                phio = 0.
                decaytime = 1000.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)
                # print 'params:', s.get_fit_params()
                # print 'errors:', s.get_fit_errors()
                plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )
                # plot2d_1.add(T_vec*1e9, values_from_fit)
                plot2d_1.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')
                # print 'pi pulse time is : '+ str(np.pi/p[2]) +' ns'

                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, imag_a)
                # guess parameters##########################################################
                background = (imag_a.max() + imag_a.min() )/2.
                osc_amp = (imag_a.max() - imag_a.min() )/2.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)

                plot2d_4.replace_inline_data_y2(T_vec*1e9, imag_a, values_from_fit )
                plot2d_4.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')

                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, real_a)
                # guess parameters##########################################################
                background = (real_a.max() + real_a.min() )/2.
                osc_amp = (real_a.max() - real_a.min() )/2.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)

                plot2d_3.replace_inline_data_y2(T_vec*1e9, real_a, values_from_fit )
                plot2d_3.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')
                print '3'
            if FIT:
                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, amplitude)
                # guess parameters##########################################################
                phio = 0.
                decaytime = 1000.

                #background = (amplitude.max() + amplitude.min() )/2.
                background = sum(amplitude)/len(amplitude)
                osc_amp = (amplitude.max() - amplitude.min() )/2.

                #background_pha = (phase.max() + phase.min() )/2.
                background_pha = sum(phase)/len(phase)
                osc_pha = (phase.max() - phase.min() )/2.

                # fitting ##################################################################
                if nom_expect != 0:
                    nb_expected_oscillation = nom_expect
                    pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
                else:
                #!V 180824 automaticaly counts Rabi-pulses
                    nom_cross = 0
                    for i in np.arange(3, len(amplitude)):
                        #s[i] positive sign or not
                        s0 = amplitude[i-3] > background
                        s1 = amplitude[i-2] > background
                        s2 = amplitude[i-1] > background
                        s3 = amplitude[i]   > background
                        if (s0==s1 and s2==s3 and s1!=s3):
                            nom_cross = nom_cross +1
                    nb_expected_oscillation = nom_cross/2
                    #print nb_expected_oscillation
                    pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation

                if Ajust_window: #now switched off
                    if cycle_counter > 3:
                        if nb_expected_oscillation > 5:
                            coff = nb_expected_oscillation/5
                            window0 = window/coff
                            tpi = get_rabi_pi(parameters, power2, GausSigma=GausSigma, window = window0)
                            return tpi

                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)
                # print 'params:', s.get_fit_params()
                # print 'errors:', s.get_fit_errors()
                plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )
                # plot2d_1.add(T_vec*1e9, values_from_fit)
                plot2d_1.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')
                tpi = np.pi/p[2]
                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, imag_a)
                # guess parameters##########################################################
                background = (imag_a.max() + imag_a.min() )/2.
                osc_amp = (imag_a.max() - imag_a.min() )/2.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)

                plot2d_4.replace_inline_data_y2(T_vec*1e9, imag_a, values_from_fit )
                plot2d_4.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')

                s = fit.ExponentialDecaySine()
                s.set_data(T_vec*1e9, real_a)
                # guess parameters##########################################################
                background = (real_a.max() + real_a.min() )/2.
                osc_amp = (real_a.max() - real_a.min() )/2.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)

                plot2d_3.replace_inline_data_y2(T_vec*1e9, real_a, values_from_fit )
                plot2d_3.set_plottitle('T_pi= '+ str(np.pi/p[2])+' ns'+', T_rabi= '+str(p[4])+' ns')
        print '3' #!V
        Tabor.set_trigger_source('EVEN')
        if Pulsing_instrument.get_board_flag():
            Pulsing_instrument.measurement_close(transfert_info=False)

    except:
        print 'WARNING: (rabi) An error ocurred:'
        e = sys.exc_info()[1]
        print e.args[0]
    finally:
        if Pulsing_instrument.get_board_flag():
            Pulsing_instrument.measurement_close(transfert_info=False)

        data_measurement.add_data_point(T_vec*1e9, amplitude, phase, real_a, imag_a)
        data_measurement.close_file()
        if COMMENTS:
            print Pulsing_instrument.measurement_close(transfert_info=True)
        else:
            Pulsing_instrument.measurement_close(transfert_info=True)

        Tabor.set_trigger_source('EVEN')


    if FIT:
        data_fit.add_data_point( s.get_fit_params(), s.get_fit_errors())
        data_fit.close_file()
        plot2d_1.add(T_vec*1e9, values_from_fit)



    plot2d_1.save_png()
    plot2d_2.save_png()
    plot2d_3.save_png()
    plot2d_4.save_png()

    # this is empirical...
    if COMMENTS:
        print 'distance=', round(np.sqrt((real_a[9]-real_a[0])**2 + (imag_a[9]-imag_a[0])**2)*1e3, 2), T_vec[9]*1e9


    qt.mend()
    return tpi*1e-9

def get_ramsey_df(parameters, Tabor_loading=1, averaging = 3e3):
    if COMMENTS:
        print 'Ramsey, freq-q=', parameters.freq_q
    #returns frequency diference to f_q
    df=0
    FIT = True
    try:
        power1 = parameters.power1
        power2 = parameters.power2
        f_atom = parameters.freq_q
        f_cav = parameters.freq_read
        t_meas = parameters.t_read
        tpi = parameters.tpi
        GausSigma = parameters.nsigma
        if GausSigma !=0:
            print 'Warning!'
            print 'Sorry, but Ramsey better work with square pulses. Choose another power2 and Tpi please'
            print 'Ramsey run anyway with correct pait of power2 and tpi, but pi_o_2 pulses are not perfect!'
            # write a program of calibrating tpi_o_2 as tpi
            print 'Calibrate freq_q with square pulses and after switch to gaussian'
    except:
        print 'WARNING: (ramsey) one of parameters is not defined!'
        return None

    #t_pi_o2 = 1e-9 * math.floor(tpi/2) # we define the excitation time for a pi over 2 pulse on the qubit
    t_pi_o2 = tpi/2 #should be in sec
    if COMMENTS:
        print 't_pi_o2: ', t_pi_o2
    delta_t = 200 # not to touch for now
    acq_time = t_meas*1e9 + delta_t + 300. # not to touch for now

    # We define the time vector of the Ramsey measurement:
    tr_stop = 1e-6 # in s
    tr_step = 5e-9 # in s
    tr_start = 0e-6 # in s
    T_vec = np.arange(tr_start, tr_stop, tr_step)

    nb_sequences = len(T_vec) # we define the number of sequences to acquire for the board.

    t_wait = 0e-9
    t_rise = None
    ###########################################################
    #
    #
    #               Experiment
    #
    #
    ###########################################################

    if Tabor_loading:
        # here we write in the AWG memory if the condition Tabor_loading is True
        if COMMENTS:
            print 'Tabor loading...'

        if GausSigma < 0:
            GausSigma = 0
            print 'Warning_Fidelity: Gaussigma was <0, changed to 0'
        if GausSigma > 5:
            GausSigma = 5
            print 'Warning_Fidelity: Gaussigma was >5, changed to 5'

        Pulsing_instrument.set_trigger_time(100.)
        Pulsing_instrument.write_Ramsey_pulsessequence(t_pi_o2, tr_stop, tr_step,
                tr_start, t_meas, t_wait=t_wait, delete='all', t_rise =t_rise, nsigma=GausSigma )


    qt.mstart()

    data_measurement = qt.Data(name='Ramsey')
    data_measurement.add_coordinate('waiting time [ns]', units = 'ns')
    data_measurement.add_value('S21 ',            units = 'Volt')
    data_measurement.add_value('Phase ',            units = 'rad')
    data_measurement.add_value('Re ',            units = 'Volt')
    data_measurement.add_value('Im ',            units = 'Volt')
    data_measurement.create_file()


    plot2d_1 = qt.Plot2D(data_measurement,
                      name      = 'S21 ramsey',
                      coorddim  = 0,
                      valdim    = 1,
                      maxtraces = 2)

    plot2d_2 = qt.Plot2D(data_measurement,
                        name      = 'Phase ramsey',
                        coorddim  = 0,
                        valdim    = 2,
                        maxtraces = 2)
    #
    plot2d_3 = qt.Plot2D(data_measurement,
                      name      = 'Re ramsey',
                      coorddim  = 0,
                      valdim    = 3,
                      maxtraces = 2)

    plot2d_4 = qt.Plot2D(data_measurement,
                        name      = 'Im ramsey',
                        coorddim  = 0,
                        valdim    = 4,
                        maxtraces = 2)

    if FIT:

        data_fit = qt.Data(name='Ramsey_OSC_fit')
        data_fit.add_value('parameters ',            units = 'rad, rad, GHz*2pi, rad, ns')
        data_fit.add_value('errors ',            units = 'rad, rad, GHz*2pi, rad, ns')
        data_fit.create_file()

    try:
        Pulsing_instrument.prep_rabi(f_cav, f_atom, averaging, nb_sequences,
            power1, power2, acq_time, t_meas*1e9, delta_t)
        # Pulsing_instrument.prep_ramsey(f_cav, f_atom, averaging, nb_sequences, power1, power2)
        qt.msleep(3)
        Tabor.set_trigger_source('TIM')
        while Pulsing_instrument.get_acquisition_completed() != 100.:
            #print  Pulsing_instrument.get_acquisition_completed(), '%'

            result = Pulsing_instrument.measurement()
            ((real_a, rea0), (imag_a, ima0))= result
            real_a -= rea0
            imag_a -= ima0
            amplitude = np.sqrt(real_a**2+imag_a**2)

            complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)

            # (real_a, imag_a)= result
            #
            # amplitude = np.sqrt(real_a**2+imag_a**2)
            #
            # complexe = (real_a + 1j*imag_a )*np.exp(1j*f_cav*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            # phase = np.angle(complexe)

            qt.msleep(0.1)

            # plot2d_1.replace_inline_data(T_vec*1e9, amplitude)
            plot2d_2.replace_inline_data(T_vec*1e9, phase)
            # plot2d_3.replace_inline_data(T_vec*1e9, real_a)
            # plot2d_4.replace_inline_data(T_vec*1e9, imag_a)
            if FIT:
                s = fit.ExponentialDecaySine()
                signal = real_a
                s.set_data(T_vec*1e9, signal)
                # guess parameters##########################################################
                background = (signal.max() + signal.min() )/2.
                osc_amp = (signal.max() - signal.min() )/2.
                nb_expected_oscillation = 4.
                pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
                phio = 0.
                decaytime = 1000.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)
                plot2d_3.replace_inline_data_y2(T_vec*1e9, signal, values_from_fit )
                plot2d_3.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')

                s = fit.ExponentialDecaySine()
                signal = imag_a
                s.set_data(T_vec*1e9, signal)
                # guess parameters##########################################################
                background = (signal.max() + signal.min() )/2.
                osc_amp = (signal.max() - signal.min() )/2.
                nb_expected_oscillation = 10.
                pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
                phio = 0.
                decaytime = 1000.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)
                plot2d_4.replace_inline_data_y2(T_vec*1e9, signal, values_from_fit )
                plot2d_4.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')

                s = fit.ExponentialDecaySine()
                signal = amplitude
                s.set_data(T_vec*1e9, signal)
                # guess parameters##########################################################
                background = (signal.max() + signal.min() )/2.
                osc_amp = (signal.max() - signal.min() )/2.
                nb_expected_oscillation = 10.
                pulsation = 2*np.pi/(tr_stop*1e9)*nb_expected_oscillation
                # print pulsation/2/np.pi
                phio = 0.
                decaytime = 1000.
                p0 = [background, osc_amp, pulsation, phio, decaytime]
                # fitting ##################################################################
                p = s.fit(p0)
                values_from_fit = s.func(p)
                # print 'params:', s.get_fit_params()
                # print 'errors:', s.get_fit_errors()

                plot2d_1.replace_inline_data_y2(T_vec*1e9, amplitude, values_from_fit )

                # plot2d_1.add(T_vec*1e9, values_from_fit)
                plot2d_1.set_plottitle('T2= '+ str(p[4])+' ns'+', Freq= '+str(p[2]/2/np.pi*1e3)+' MHz')
                #df = p[2]/2/np.pi*1e3
                #print df

        Tabor.set_trigger_source('EVEN')
        if Pulsing_instrument.get_board_flag():
            Pulsing_instrument.measurement_close(transfert_info=False)

        #df = p[2]/2/np.pi*1e3
        df = p[2]/2/np.pi
        #print df

    finally:
        if Pulsing_instrument.get_board_flag():
            Pulsing_instrument.measurement_close(transfert_info=False)

        data_measurement.add_data_point(T_vec*1e9, amplitude, phase, real_a, imag_a)
        data_measurement.close_file()

        #print Pulsing_instrument.measurement_close(transfert_info=True)
        Pulsing_instrument.measurement_close(transfert_info=True)
        # smb_atom.set_freqsweep('OFF')
        # smb_atom.set_gui_update('ON')
        Tabor.set_trigger_source('EVEN')

        if FIT:
            plot2d_1.add(T_vec*1e9, values_from_fit)
            data_fit.add_data_point( s.get_fit_params(), s.get_fit_errors())
            data_fit.close_file()


    plot2d_1.save_png()
    plot2d_2.save_png()
    qt.mend()
    if COMMENTS:
        print 'delta f = ', "%.8f" %df, 'GHz'
    df_mhz = df*1e3
    if COMMENTS:
        print 'delta_f = ', "%.4f" %df_mhz, 'MHz'
    return df

def get_dblob(f_min, f_max, parameters, averaging = 1e3, nop = 100, Tabor_loading = 1):
    '''
    Makes sweep s21/ph vs readout frequences, returns optimal freq_read
    transmission(f_min, f_max, averaging = 1e3, nop = 100, Tabor_loading = 1)
    - now just maximum of d-blobs (how to fit it?)
    [just repeat in small spot with more av]
    '''
    try:
        power1 = parameters.power1
        power2 = parameters.power2
        freq_q = parameters.freq_q
        freq_cav = parameters.freq_read
        tpi = 1e-9 * parameters.tpi
        t_read = parameters.t_read
    except:
        print 'W@RNING: (dblob) one of parameters is not defined!'
        return None
    #start
    if COMMENTS:
        print 'Dblob, f_read=', parameters.freq_read
    #Step of the Sweep [GHz]
    f_step = (abs(f_max-f_min)/nop)
    if COMMENTS:
        print('fstep=', f_step)

    freq_vec = np.arange(f_min, f_max + f_step, f_step)
    if len(freq_vec) %2 !=0:
            freq_vec = np.arange(f_min, f_max , f_step)

    #tpi = 40e-9
    #t_read = 500e-9

    delta_t = 0.2e-6
    acq_time =  t_read*1e9 + delta_t*1e9 + 300.
    t_rise = None
    # 10e-9
    tau = 50.
    ###########################################################
    #
    #
    #               Experiment
    #
    #
    ###########################################################
    if Tabor_loading:
        if COMMENTS:
            print 'Tabor loading...'
        Pulsing_instrument.set_trigger_time(100.)
        Pulsing_instrument.write_twotone_pulsessequence_withpi(temp_1=t_read,
            t1_start= tpi + 0.1e-6, temp_2=tpi , m1_start= tpi,  delete = 'all', t_rise=t_rise)


    # Pulsing_instrument.write_twotone_pulsessequence( 500e-9, 100e-9 + tpi, tpi, delete = 'all')
    qt.mstart()

    data_measurement = qt.Data(name='Cavity_shift')
    data_measurement.add_coordinate('R.O. frequency [GHz]', units = 'GHz')
    data_measurement.add_value('S21 ',            units = 'Volt')
    data_measurement.add_value('Phase ',            units = 'rad')
    data_measurement.add_value('Re ',            units = 'Volt')
    data_measurement.add_value('Im ',            units = 'Volt')
    data_measurement.add_value('S21 pi',            units = 'Volt')
    data_measurement.add_value('Phase pi',            units = 'rad')
    data_measurement.add_value('Re pi',            units = 'Volt')
    data_measurement.add_value('Im pi',            units = 'Volt')
    data_measurement.add_value('D blobs',            units = 'Volt')

    data_measurement.create_file()

    # plot2d_1 = qt.Plot2D(data_measurement,
    #                   name      = 'S21 ',
    #                   coorddim  = 0,
    #                   valdim    = 1,
    #                   maxtraces = 2)
    # #
    # plot2d_2 = qt.Plot2D(data_measurement,
    #                     name      = 'Phase ',
    #                     coorddim  = 0,
    #                     valdim    = 2,
    #                     maxtraces = 2)
    # #
    # plot2d_3 = qt.Plot2D(data_measurement,
    #                   name      = 'Re ',
    #                   coorddim  = 0,
    #                   valdim    = 3,
    #                   maxtraces = 2)
    #
    # plot2d_4 = qt.Plot2D(data_measurement,
    #                     name      = 'Im ',
    #                     coorddim  = 0,
    #                     valdim    = 4,
    #                     maxtraces = 2)

    plot2d_5 = qt.Plot2D(data_measurement,
                        name      = 'd blobs ',
                        coorddim  = 0,
                        valdim    = 9,
                        maxtraces = 2)
    # With Pi-Pulse:
    board_flag = None
    #d_dict = {}
    try:
        # Pulsing_instrument.prep_conditional_transmission(freq_vec, averaging,
        #             power1, f_cw=freq_q, power2=power2, acq_time=acq_time, pulse_time=t_read*1e9, delta_t=delta_t )
        Pulsing_instrument.prep_conditional_transmission(freq_vec, averaging,
                    power1, f_cw=freq_q, power2=power2, acq_time=acq_time, pulse_time=t_read*1e9, delta_t=delta_t, tau=tau )

        qt.msleep(2)
        # smb_atom.set_freqsweep('OFF')
        smb_cavity.restartsweep()
        qt.msleep(1)

        board_flag = True
        Tabor.set_trigger_source('TIM')
        while ats9360.get_completed_acquisition() != 100.:
            #print  ats9360.get_completed_acquisition(), '%'
            result = ats9360.measurement()
            # (real, imag)= result
            ((real, rea0), (imag,ima0))= result
            real = real - np.mean(rea0)
            imag = imag - np.mean(ima0)

            real = np.reshape(real, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) )
            imag = np.reshape(imag, (len(freq_vec), Pulsing_instrument.get_pulsenumber_averaging()) )
            real_a = real[:, :Pulsing_instrument.get_pulsenumber_averaging()/2]
            imag_a = imag[:, :Pulsing_instrument.get_pulsenumber_averaging()/2]
            real_a = np.mean(real_a, axis = 1)
            imag_a = np.mean(imag_a, axis = 1)
            real_a_pi = real[:, Pulsing_instrument.get_pulsenumber_averaging()/2:]
            imag_a_pi = imag[:, Pulsing_instrument.get_pulsenumber_averaging()/2:]
            real_a_pi = np.mean(real_a_pi, axis = 1)
            imag_a_pi = np.mean(imag_a_pi, axis = 1)

            d_blobs = np.sqrt((real_a-real_a_pi)**2+(imag_a-imag_a_pi)**2)

            amplitude = np.sqrt(real_a**2+imag_a**2)
            complexe = (real_a + 1j*imag_a )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase = np.angle(complexe)

            amplitude_pi = np.sqrt(real_a_pi**2+imag_a_pi**2)
            complexe_pi = (real_a_pi + 1j*imag_a_pi )*np.exp(1j*freq_vec*Pulsing_instrument.get_electrical_phase_delay()*2*np.pi)
            phase_pi = np.angle(complexe_pi)
            qt.msleep(0.1)
            #
        #    plot2d_1.replace_inline_data_y2(freq_vec, amplitude, amplitude_pi)
        #    plot2d_2.replace_inline_data_y2(freq_vec, phase, phase_pi)
        #    plot2d_3.replace_inline_data_y2(freq_vec, real_a, real_a_pi)
        #    plot2d_4.replace_inline_data_y2(freq_vec, imag_a, imag_a_pi)
            plot2d_5.replace_inline_data(freq_vec, d_blobs)

        Tabor.set_trigger_source('EVEN')
        ats9360.measurement_close(transfert_info=False)
        board_flag = False
    finally:
        if board_flag:
            ats9360.measurement_close(transfert_info=False)

        data_measurement.add_data_point(freq_vec,amplitude,phase, real_a, imag_a, amplitude_pi, phase_pi, real_a_pi, imag_a_pi, d_blobs)
    #    plot2d_1.add(freq_vec, amplitude_pi)
    #    plot2d_2.add(freq_vec, phase_pi)
    #    plot2d_3.add(freq_vec, real_a_pi)
    #    plot2d_4.add(freq_vec, imag_a_pi)

        #print ats9360.measurement_close(transfert_info=True)
        # smb_cavity.set_freqsweep('OFF')
        # smb_cavity.set_gui_update('ON')
        Tabor.set_trigger_source('EVEN')

    #return max[d_blobs]
    data_measurement.close_file()
    #plot2d_1.save_png()
    #plot2d_2.save_png()
    # plot2d_3.save_png()
    # plot2d_4.save_png()
    plot2d_5.save_png()

    qt.mend()
    # Searching for maximum:
    best_dblob = np.max(d_blobs)
    best_freq = freq_vec[np.argmax(d_blobs)]
    result = {'best_freq':best_freq,'best_dblob':best_dblob}
    return result

def get_fidelity(parameters, Tabor_loading = 1, counts = 10000, average = 12):
    '''
    Starts SingleShot readout with postselection.
    Returns dictionary {Fidelities(F, Fleft, Fright), Fidelities_postselected(F, Fg, Fe), Fidelities_gaussian(F, 1-eg, 1-ge), errors_without_owerlap(e,g)}
    ## keep this number of average even and below 20. no need for too big data set for histograms.
    '''
    FIT = False
    Fidelity = 1

    try:
        power1 = parameters.power1
        power2 = parameters.power2
        freq_q = parameters.freq_q
        freq_cav = parameters.freq_read
        tpi = parameters.tpi
        t_read = parameters.t_read
        GausSigma = parameters.nsigma
    except:
        print 'WARNING: (Fidelity) one of parameters is not defined!'
        return None
    #start
    print 'SingleShot...'

    t1_1 = t_read #?for first readout in measurment (g)
    t1_2 = t_read #?for second readout in measurment (g or e)

    t_rise = None
    # t_rise = 10e-9
    tau = 50
    # 65.
    # 30.
    delta_t = 250.

    # acq_time = 2600.
    t1_1start = 0.1e-6
    t_between = 200e-9
    t2_start = t1_1start + t1_1 + t_between
    t1_2start = t2_start + tpi

    t_protect = 0.

    acq_time = (t1_2start + t1_2)*1e9 + delta_t + 300.

    delta_m1_start = 0.1e-6
    dmstart = 0.1e-6
    threshold = 0e-3 # here signal is following imag part

    mw = 2
    Pulsing_instrument.set_routing_awg({'secondtone_channel':4})

    if Tabor_loading:
        if COMMENTS:
            print 'Tabor loading...'
        Pulsing_instrument.set_trigger_time(200.)
        if GausSigma < 0:
            GausSigma = 0
            print 'Warning_Fidelity: Gaussigma was <0, changed to 0'
        if GausSigma > 5:
            GausSigma = 5
            print 'Warning_Fidelity: Gaussigma was >5, changed to 5'
        #gaussian pulse here                                                                #    |-\___
        #diff_PdB=-0 - means that readout pulse is simple, but if you change it pulse become ___|      \_
        Pulsing_instrument.write_IQpi_several_RO_bifurcationshape(tpi, t1_1, t1_2, duty_cycle=0.05, diff_PdB=-0, t2_start=t2_start ,
                t1_1start=t1_1start, t1_2start=t1_2start, t_rise= t_rise, delta_m1_start=delta_m1_start, delete=True, nsigma=GausSigma)

    qt.mstart()

    data_measurement = qt.Data(name='IQ_postselection_t1='+str(t_read*1e9))
    data_measurement.add_coordinate('Real ',         units = 'volt')
    data_measurement.add_coordinate('Real-pi',         units = 'volt')
    data_measurement.add_value('Imag ',         units = 'Volt')
    data_measurement.add_value('Imag-pi ',         units = 'Volt')
    data_measurement.add_coordinate('Real preselect ',         units = 'volt')
    data_measurement.add_coordinate('Real-pi preselect',         units = 'volt')
    data_measurement.add_value('Imag preselect ',         units = 'Volt')
    data_measurement.add_value('Imag-pi preselect',         units = 'Volt')

    data_measurement.create_file()

    data_measurement_hist= qt.Data(name='plot_hist')
    data_measurement_hist.add_coordinate('Real ',         units = 'volt')
    data_measurement_hist.add_coordinate('Imag ',         units = 'volt')
    data_measurement_hist.add_coordinate('real before ',         units = 'volt')
    data_measurement_hist.add_coordinate('Imag  before',         units = 'rad')
    # data_measurement_hist.add_coordinate('real postselected ',         units = 'volt')
    # data_measurement_hist.add_coordinate('Imag  postselected',         units = 'rad')

    data_measurement_hist.add_value('counts_re')
    data_measurement_hist.add_value('counts_im')
    data_measurement_hist.add_value('counts re before')
    data_measurement_hist.add_value('counts_im before')
    # data_measurement_hist.add_value('counts re postselected')
    # data_measurement_hist.add_value('counts_im postselected')

    data_measurement_hist.create_file()
    plothist_re = qt.Plot2D(data_measurement_hist,
                            name      = 'Hist Re postsel',
                            coorddim  = 0,
                            valdim    = 4)
    plothist_re.set_style('steps')
    plothist_re.set_ylog(True)
    plothist_im = qt.Plot2D(data_measurement_hist,
                            name      = 'Hist Im postsel',
                            coorddim  = 1,
                            valdim    = 5 )
    plothist_im.set_style('steps')
    plothist_im.set_ylog(True)

    # Pi pulse
    board_flag = None
    try:
        real_a = []
        imag_a = []
        real_b = []
        imag_b = []

        Pulsing_instrument.prep_IQ_2_sevRO(counts, average, freq_cav, power1,
                freq_q, power2, acq_time, t1_1*1e9 - t_protect, 0.,
                t1_2*1e9 - t_protect, t1_2start*1e9 -t1_1start*1e9, delta_t, tau=tau)
        # Pulsing_instrument.prep_IQ_2_sevRO(counts, average, freq_read, power1,
        #         freq_q, power2, acq_time, t1_1*1e9 - t_protect, t1_1start*1e9  -100,
        #         t1_2*1e9 - t_protect, t1_2start*1e9 -100, delta_t, tau=50.)
        qt.msleep(2)
        board_flag = True
        Tabor.set_trigger_source('TIM')
        while ats9360.get_completed_acquisition() != 100.:
            # print 'percentage done:' + str(ats9360.get_completed_acquisition())+'%'

            result = ats9360.measurement()

            ((real, real2, re0), (imag, imag2, im0)) = result
            # print np.shape(real)
            real_a = np.append(real_a, real - np.mean(re0))
            imag_a = np.append(imag_a, imag - np.mean(im0))
            real_b = np.append(real_b, real2 - np.mean(re0))
            imag_b = np.append(imag_b, imag2 - np.mean(im0))
            # real_a = np.append(real_a, real-rea0)
            # imag_a = np.append(imag_a, imag-ima0)
            #print np.shape(real_a) #!V


            # plot2d_2.replace_inline_data(real_a, imag_a)

            qt.msleep(0.1)

        Tabor.set_trigger_source('EVEN')
        ats9360.measurement_close(transfert_info=False)
        board_flag = False
    finally:
        if board_flag:
            ats9360.measurement_close(transfert_info=False)

        ####
        re_1 = real_b[::2]
        re_2 = real_b[1::2]
        im_1 = imag_b[::2]
        im_2 = imag_b[1::2]
        h_re = np.histogram(re_1, bins = 100, density=0)
        h_im = np.histogram(im_1, bins = 100, density=0)
        h_re2 = np.histogram(re_2, bins = 100, density=0)
        h_im2 = np.histogram(im_2, bins = 100, density=0)
        re_hist_1 = (h_re[1][:len(h_re[0])], h_re[0])
        im_hist_1 = (h_im[1][:len(h_im[0])], h_im[0])
        re_hist_2 = (h_re2[1][:len(h_re2[0])], h_re2[0])
        im_hist_2 = (h_im2[1][:len(h_im2[0])], h_im2[0])

        #####
        re_1_before = real_a[::2]
        re_2_before = real_a[1::2]
        im_1_before = imag_a[::2]
        im_2_before = imag_a[1::2]
        h_re_before = np.histogram(re_1_before, bins = 100, density=0)
        h_im_before = np.histogram(im_1_before, bins = 100, density=0)
        h_re2_before = np.histogram(re_2_before, bins = 100, density=0)
        h_im2_before = np.histogram(im_2_before, bins = 100, density=0)
        re_hist_1_before = (h_re_before[1][:len(h_re_before[0])], h_re_before[0])
        im_hist_1_before = (h_im_before[1][:len(h_im_before[0])], h_im_before[0])
        re_hist_2_before = (h_re2_before[1][:len(h_re2_before[0])], h_re2_before[0])
        im_hist_2_before = (h_im2_before[1][:len(h_im2_before[0])], h_im2_before[0])

        #####

        #print ats9360.measurement_close(transfert_info=True) #!V (no printing)
        ats9360.measurement_close(transfert_info=False)
        Tabor.set_trigger_source('EVEN')

    data_measurement.add_data_point(re_1, re_2, im_1, im_2, re_1_before , re_2_before , im_1_before , im_2_before )
    # if len(re_1_new)>0 and len(re_2_new)>0:
    #     data_measurement_post1.add_data_point(re_1_new,  im_1_new)
    #     data_measurement_post2.add_data_point(re_2_new,  im_2_new)

    # plot2d_1.add(re_2, im_2)
    data_measurement_hist.add_data_point(re_hist_1[0],  im_hist_1[0], re_hist_1_before[0], im_hist_1_before[0],# re_hist_1_new[0], im_hist_1_new[0],
                    re_hist_1[1], im_hist_1[1],re_hist_1_before[1], im_hist_1_before[1])#, re_hist_1_new[1], im_hist_1_new[1])
    data_measurement_hist.new_block()
    data_measurement_hist.add_data_point(re_hist_2[0],  im_hist_2[0], re_hist_2_before[0], im_hist_2_before[0], #re_hist_2_new[0], im_hist_2_new[0],
                        re_hist_2[1], im_hist_2[1], re_hist_2_before[1], im_hist_2_before[1])#, re_hist_2_new[1], im_hist_2_new[1])


    # plothist_re.add(re_hist_1_new[0], re_hist_1_new[1])
    # plothist_re.add(re_hist_2_new[0], re_hist_2_new[1])
    plothist_re.add(re_hist_1_before[0], re_hist_1_before[1])
    plothist_re.add(re_hist_2_before[0], re_hist_2_before[1])

    # plothist_im.add(im_hist_1_new[0], im_hist_1_new[1])
    # plothist_im.add(im_hist_2_new[0], im_hist_2_new[1])

    plothist_im.add(im_hist_1_before[0], im_hist_1_before[1])
    plothist_im.add(im_hist_2_before[0], im_hist_2_before[1])
    # plothist_re.set_ylog(1)
    if FIT:
        data_fit= qt.Data(name='Spectro_fit')
        data_fit.add_value('parameters ',            units = 'none, none, GHz, Volt')
        data_fit.add_value('errors ',            units = 'none, none, GHz, Volt')
        data_fit.create_file()
        s = fit.DoubleGaussian()

        # fitting Re without pi:
        s.set_data(re_hist_1[0],  re_hist_1[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position1 = np.mean(re_hist_1[0])
        position2 = np.mean(re_hist_2[0])
        fwhm = np.std(re_hist_1[0])
        p = [background, area, position1, fwhm, area, position2, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        plothist_re.add(re_hist_1[0], values_from_fit,'-', linewidth = 1.5)
        data_fit.new_block()
        # fitting Im without pi:
        s = fit.DoubleGaussian()
        s.set_data(im_hist_1[0],  im_hist_1[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position = np.mean(im_hist_1[0])
        fwhm = np.std(im_hist_1[0])
        position2 = np.mean(im_hist_2[0])
        p = [background, area, position, fwhm, area, position2, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        plothist_im.add(im_hist_1[0], values_from_fit, '-', linewidth = 1.5)
        data_fit.new_block()

        # fitting amp without pi:
        s = fit.Gaussian()
        s.set_data(amp_hist_1[0],  amp_hist_1[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position = np.mean(amp_hist_1[0])
        fwhm = np.std(amp_hist_1[0])
        p = [background, area, position, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        # plothist_amp.add(amp_hist_1[0], values_from_fit, '-', linewidth = 1.5)
        data_fit.new_block()
        # fitting phase without pi:
        s.set_data(phase_hist_1[0],  phase_hist_1[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position = np.mean(phase_hist_1[0])
        fwhm = np.std(phase_hist_1[0])
        p = [background, area, position, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'phase params:', s.get_fit_params()
        # print 'phase errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        # plothist_phase.add(phase_hist_1[0], values_from_fit, '-', linewidth = 1.5)
        data_fit.new_block()
        # phase_mean = np.append(phase_mean, p[2] )


        # fitting Re with pi:
        s = fit.DoubleGaussian()
        s.set_data(re_hist_2[0],  re_hist_2[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position1 = np.mean(re_hist_1[0])
        position2 = np.mean(re_hist_2[0])
        fwhm = np.std(re_hist_2[0])
        p = [background, area, position1, fwhm, area, position2, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        plothist_re.add(re_hist_2[0], values_from_fit, '-', linewidth = 1.5)
        data_fit.new_block()
        # fitting Im with pi:
        s = fit.DoubleGaussian()
        s.set_data(im_hist_2[0],  im_hist_2[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position = np.mean(im_hist_2[0])
        fwhm = np.std(im_hist_2[0])
        position2 = np.mean(im_hist_1[0])
        p = [background, area, position, fwhm/2., 0.3, position2, fwhm/2.]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        plothist_im.add(im_hist_2[0], values_from_fit, '-', linewidth = 1.5)

        # fitting amp with pi:
        s = fit.Gaussian()
        s.set_data(amp_hist_2[0],  amp_hist_2[1])
        # guess parameters##########################################################
        background = 0.
        area = 1.
        position = np.mean(amp_hist_2[0])
        fwhm = np.std(amp_hist_2[0])
        p = [background, area, position, fwhm]
        # fitting ##################################################################
        p = s.fit(p, fixed=[0, 1])
        values_from_fit = s.func(p)
        # print 'params:', s.get_fit_params()
        # print 'errors:', s.get_fit_errors()
        data_fit.add_data_point(s.get_fit_params(),s.get_fit_errors())
        # plothist_amp.add(amp_hist_2[0], values_from_fit, '-', linewidth = 1.5)

        data_fit.close_file()
    #plot2d_2.set_xrange(-1e-2, 1e-2)
    #plot2d_2.set_yrange(-1e-2, 1e-2)

    #plot2d_1.save_png()
    #plot2d_2.save_png()
    plothist_re.save_png()
    plothist_im.save_png()
    # plothist_amp.save_png()
    # plothist_phase.save_png()

    #added !V 180717
    #plothist_H.save_png()
    #plothist_Hpi.save_png()

    data_measurement.close_file()
    # data_measurement_post1.close_file()
    # data_measurement_post2.close_file()

    data_measurement_hist.close_file()
    #print 'distance: ', np.sqrt( (np.mean(re_1)-np.mean(re_2))**2 + (np.mean(im_1)-np.mean(im_2))**2 ) #!V
    # print 'distance: ', np.sqrt( (np.mean(re_1_new)-np.mean(re_2_new))**2 + (np.mean(im_1_new)-np.mean(im_2_new))**2 )

    # print np.mean(amplitude1), np.mean(amplitude2)
    qt.mend()

    # !V 180717 add parameters
    data_measurement_hist2D = qt.Data(name='plot_hist2D'+parameters.get_string())
    data_measurement_hist2D.add_coordinate('X ',         units = 'millivolt')
    data_measurement_hist2D.add_coordinate('Y ',         units = 'millivolt')
    data_measurement_hist2D.add_coordinate('X pi ',         units = 'volt')
    data_measurement_hist2D.add_coordinate('Y pi',         units = 'rad')

    data_measurement_hist2D.add_value('counts H')
    data_measurement_hist2D.add_value('counts H pi')

    data_measurement_hist2D.create_file()

    plothist_H = qt.Plot3D(data_measurement_hist2D,
                            name      = 'H ',
                            coorddim  = (0,1),
                            valdim    = 4 )
    plothist_H.set_palette('bluewhitered')
    plothist_Hpi = qt.Plot3D(data_measurement_hist2D,
                            name      = 'H pi',
                            coorddim  = (2,3),
                            valdim    = 5)
    plothist_Hpi.set_palette('bluewhitered')

    # added by !V 180717
    plothist_H.save_png()
    plothist_Hpi.save_png()

    if Fidelity:

        Re = re_1
        Im = im_1
        Re_pi = re_2
        Im_pi = im_2
        Vmax = np.max((Re, Im, Re_pi, Im_pi))*1e3
        Vmin = np.min((Re, Im, Re_pi, Im_pi))*1e3
        H, xe, ye = np.histogram2d(Re, Im, bins=(100, 100))
        H = H.T
        X, Y = np.meshgrid(xe*1e3, ye*1e3)

        H_pi, xe_pi, ye_pi = np.histogram2d(Re_pi, Im_pi, bins=(100, 100))
        H_pi = H_pi.T
        X_pi, Y_pi = np.meshgrid(xe_pi*1e3, ye_pi*1e3)
        plothist_H.set_xrange(Vmin, Vmax)
        plothist_H.set_yrange(Vmin, Vmax)
        plothist_Hpi.set_xrange(Vmin, Vmax)
        plothist_Hpi.set_yrange(Vmin, Vmax)
        #print 'V', Vmin, Vmax #!V
        for i in np.arange(len(X[0,:])):
            if i != len(X[0,:])-1:
                data_measurement_hist2D.add_data_point(X[:-1,i], Y[:-1,i], X_pi[:-1,i], Y_pi[:-1,i], H[:,i], H_pi[:,i] )
                data_measurement_hist2D.new_block()

        data_measurement_hist2D.close_file()

        Re_post = re_1_before
        Im_post = im_1_before
        Re_pi_post = re_2_before
        Im_pi_post = im_2_before

        C1 = Re + 1j*Im
        C2 = Re_pi +1j*Im_pi
        theta = -np.angle(np.mean(C1)-np.mean(C2))
        Re = np.real(C1*np.exp(1j*theta))
        Im = np.imag(C1*np.exp(1j*theta))
        Re_pi = np.real(C2*np.exp(1j*theta))
        Im_pi = np.imag(C2*np.exp(1j*theta))

        # R = np.hstack((Re, Re_pi))
        # I = np.hstack((Im, Im_pi))


        # fig, ax = plt.subplots(1,1)
        # c = ax.pcolormesh(X, Y, H)
        # c.set_cmap('hot')
        # cc = plt.colorbar(c)
        # cc.set_cmap('hot')
        # cc.set_label('Counts')
        # plt.savefig('Hist2D_split.png',
        #             dpi=1000,
        #             bbox_inches='tight',
        #             )
        # plt.show()

        C1_post = Re_post + 1j*Im_post
        C2_post = Re_pi_post +1j*Im_pi_post

        # theta= np.angle(np.mean(C2_post)-np.mean(C1_post))

        Re_post = np.real(C1_post*np.exp(1j*theta))
        Im_post= np.imag(C1_post*np.exp(1j*theta))
        Re_pi_post = np.real(C2_post*np.exp(1j*theta))
        Im_pi_post = np.imag(C2_post*np.exp(1j*theta))

        re_hist = np.histogram(Re, bins=100,normed=1, density=1)
        re_pi_hist = np.histogram(Re_pi, bins=100, normed=1, density=1)
        im_hist = np.histogram(Im, bins=100, normed=1,density=1)
        im_pi_hist = np.histogram(Im_pi, bins=100,normed=1, density=1)

        re_hist_post = np.histogram(Re_post, bins=100,normed=1, density=1)
        re_pi_hist_post = np.histogram(Re_pi_post, bins=100, normed=1, density=1)
        im_hist_post = np.histogram(Im_post, bins=100, normed=1,density=1)
        im_pi_hist_post = np.histogram(Im_pi_post, bins=100,normed=1, density=1)

        ind1 = np.where(Re_post > threshold)
        ind2 = np.where(Re_pi_post > threshold)

        Re_postselected = np.delete(Re, ind1)
        Re_pi_postselected = np.delete(Re_pi, ind2)
        Im_postselected = np.delete(Im, ind1)
        Im_pi_postselected = np.delete(Im_pi, ind2)

        re_hist_postselected = np.histogram(Re_postselected, bins=100,normed=1, density=1)
        re_pi_hist_postselected = np.histogram(Re_pi_postselected, bins=100, normed=1, density=1)
        im_hist_postselected = np.histogram(Im_postselected, bins=100, normed=1,density=1)
        im_pi_hist_postselected = np.histogram(Im_pi_postselected, bins=100,normed=1, density=1)

        d = np.mean(Re)- np.mean(Re_pi)
        m_im = np.mean(Im)
        m_im_pi =  np.mean(Im_pi)
        m_re = np.mean(Re)
        m_re_pi =  np.mean(Re_pi)
        std_re = np.std(Re)
        std_re_pi = np.std(Re_pi)
        std_im = np.std(Im)
        std_im_pi = np.std(Im_pi)
        s_re_g = fit.Gaussian()
        s_re_g.set_data(re_hist[1][:-1],re_hist[0])
        p0 = [0, 60e3, m_re, std_re]
        p = s_re_g.fit(p0, fixed=[0])
        re_g_th = s_re_g.func(p)

        s_re_pi_g = fit.Gaussian()
        s_re_pi_g.set_data(re_pi_hist[1][:-1],re_pi_hist[0])
        p0 = [0, 60e3, m_re_pi, std_re_pi]
        p_pi = s_re_pi_g.fit(p0, fixed=[0])
        re_pi_g_th = s_re_pi_g.func(p_pi)

        s_re_dg = fit.DoubleGaussian()
        s_re_dg.set_data(re_hist[1][:-1],re_hist[0])
        p0 = np.hstack((p, [0, p_pi[2], p_pi[3]]))

        p_dg = s_re_dg.fit(p0, fixed=[0])
        re_dg_th = s_re_dg.func(p_dg)
        #print 'heights:', p_dg[1] / p_dg[3] / np.sqrt(np.pi/2), p_dg[4] / p_dg[6] / np.sqrt(np.pi/2) #!V

        s_re_pi_dg = fit.DoubleGaussian()
        s_re_pi_dg.set_data(re_pi_hist[1][:-1],re_pi_hist[0])
        p0 = np.hstack((p_pi, [0, p[2], p[3]]))

        p_pi_dg = s_re_pi_dg.fit(p0, fixed=[0])
        re_pi_dg_th = s_re_pi_dg.func(p_pi_dg)
        #print 'heights pi:', p_pi_dg[1] / p_pi_dg[3] / np.sqrt(np.pi/2), p_pi_dg[4] / p_pi_dg[6] / np.sqrt(np.pi/2) #!V
        s = fit.DoubleGaussian()
        # !V Warning error float to integer! (np.linspace)!!!
        #print 'vec', np.linspace(np.min( (Re, Re_pi, Im, Im_pi)), np.max((Re, Re_pi, Im, Im_pi)), 1e3 )
        vec = np.linspace(np.min( (Re, Re_pi, Im, Im_pi)), np.max((Re, Re_pi, Im, Im_pi)), 1e3 )

        s.set_data(vec, vec)

        re_gfdg = s.func( p_dg)# [:4])
        re_pi_gfdg = s.func(p_pi_dg)#[:4])
        ind = np.argwhere(np.diff(np.sign(re_gfdg - re_pi_gfdg)) != 0)

        s_re_dg_postselected = fit.DoubleGaussian()
        s_re_dg_postselected.set_data(re_hist_postselected[1][:-1],re_hist_postselected[0])
        p0 = np.hstack((p, [0, p_pi[2], p_pi[3]]))
        p_dg = s_re_dg_postselected.fit(p0, fixed=[0])
        re_dg_th_postselected = s_re_dg_postselected.func(p_dg)
        #print 'heights postselected:', p_dg[1] / p_dg[3] / np.sqrt(np.pi/2), p_dg[4] / p_dg[6] / np.sqrt(np.pi/2) #!V
        print 'ok3'
        s_re_pi_dg_postselected = fit.DoubleGaussian()
        s_re_pi_dg_postselected.set_data(re_pi_hist_postselected[1][:-1],re_pi_hist_postselected[0])
        p0 = np.hstack((p_pi, [0, p[2], p[3]]))
        p_pi_dg = s_re_pi_dg_postselected.fit(p0, fixed=[0])
        re_pi_dg_th_postselected = s_re_pi_dg_postselected.func(p_pi_dg)
        #print 'heights pi postselected:', p_pi_dg[1] / p_pi_dg[3] / np.sqrt(np.pi/2), p_pi_dg[4] / p_pi_dg[6] / np.sqrt(np.pi/2) #!V

        s = fit.DoubleGaussian()
        vec = np.linspace(np.min( (Re, Re_pi, Im, Im_pi)), np.max((Re, Re_pi, Im, Im_pi)), 1e3 )
        s.set_data(vec, vec)

        re_gfdg_postselected = s.func( p_dg)# [:4])
        re_pi_gfdg_postselected = s.func(p_pi_dg)#[:4])
        print 'ok-f'
        s = fit.Gaussian()
        s.set_data(vec, vec)
        R = s.func(p_dg[:4])
        R_pi = s.func(p_pi_dg[:4])
        ind = np.argwhere(np.diff(np.sign(R - R_pi)) != 0)
        if len(ind) == 1:
            threshold = vec[ind]
            val_threshold = re_gfdg[ind]
            #print 'vec',  vec[ind] #!V
        elif len(ind)>1:
            #print 'vec',  vec[ind] #!V
            threshold0 = vec[ind[0]]
            threshold1 = vec[ind[1]]
            if abs(threshold - threshold0)<0.5e-3:
                threshold = threshold0
                val_threshold = re_gfdg[ind[0]]
            elif abs(threshold - threshold1)<0.5e-3:
                threshold = threshold1
                val_threshold = re_gfdg[ind[1]]
        # plt.plot(re_hist[1][:-1],re_hist[0])
        # plt.plot(re_pi_hist[1][:-1],re_pi_hist[0])
        # plt.yscale('log')
        # plt.grid()
        # plt.show()
        # print 'automatic threshold is:', threshold
        # threshold = input("Please, enter threshold value: ")
        # threshold = 2.3e-3
        # threshold = 1.8e-3
        # print 'threshold', threshold
        q_b = 0.
        q_a = 0.
        q_b_G = 0
        q_a_G = 0

        q_b_post = 0.
        q_a_post = 0.
        for j in np.arange(len(vec)):
            if vec[j]<threshold:
                q_b_G += R[j]
            else:
                q_a_G += R[j]
        PG_eg = q_b_G/(q_a_G+q_b_G)
        #print 'gaussian P(e_g): ', PG_eg, q_b_G, q_a_G #!V

        for j in np.arange(len(re_hist[1][:-1])):
            val = re_hist[1][j]
            if  val < threshold:
                q_b += re_hist[0][j]
            else:
                q_a += re_hist[0][j]
        P_e_g = float( q_b)/(q_a+q_b)
        #print 'P_e_g: ', P_e_g, q_a, q_b #!V
        for j in np.arange(len(re_hist_postselected[1][:-1])):
            val = re_hist_postselected[1][j]
            if  val < threshold:
                q_b_post += re_hist_postselected[0][j]
            else:
                q_a_post += re_hist_postselected[0][j]
        P_e_g_post = float( q_b_post)/(q_a_post+q_b_post)
        #print 'P_e_g postselected: ', P_e_g_post, q_a_post, q_b_post #!V
        q_b_pi_G = 0
        q_a_pi_G = 0
        for j in np.arange(len(vec)):
            if vec[j]<threshold:
                q_b_pi_G += R_pi[j]
            else:
                q_a_pi_G += R_pi[j]

        PG_ge = q_a_pi_G/(q_a_pi_G+q_b_pi_G)
        #print 'gaussian P(g_e): ', PG_ge, q_b_pi_G, q_a_pi_G #!V

        q_b_pi = 0.
        q_a_pi = 0.
        q_b_pi_post = 0.
        q_a_pi_post = 0.
        for j in np.arange(len(re_pi_hist[1][:-1])):
            val = re_pi_hist[1][j]
            if  val < threshold:
                q_b_pi += re_pi_hist[0][j]
            else:
                q_a_pi += re_pi_hist[0][j]
        P_g_e = float(q_a_pi)/(q_a_pi+q_b_pi)
        #print 'P_g_e: ', P_g_e, q_a_pi, q_b_pi #!V

        for j in np.arange(len(re_pi_hist_postselected[1][:-1])):
            val = re_pi_hist_postselected[1][j]
            if  val < threshold:
                q_b_pi_post += re_pi_hist_postselected[0][j]
            else:
                q_a_pi_post += re_pi_hist_postselected[0][j]

        P_g_e_post = float(q_a_pi_post)/(q_a_pi_post+q_b_pi_post)
        #print 'P_g_e postselected: ', P_g_e_post, q_a_pi_post, q_b_pi_post #!V

        F_RO = 1. - P_e_g - P_g_e
        F_g = 1. - P_e_g
        F_e = 1. - P_g_e
        F_RO_G = 1. - PG_eg - PG_ge
        F_g_G = 1. - PG_eg
        F_e_G = 1. - PG_ge
        F_RO_post = 1. - P_e_g_post - P_g_e_post
        F_e_P = 1- P_g_e_post
        F_g_P = 1-P_e_g_post

        #!V (no printing)
        # print '--------------'
        # print 'Fidelities:'
        # print '--------------'
        # print 'readout fidelities= ', F_RO, 'Fleft= ', F_g, 'Fright= ', F_e
        # print 'fidelities postselected: ', F_RO_post, F_g_P, F_e_P
        # print 'fidelities Gaussian: ', 1.-PG_eg-PG_ge, 1.-PG_eg, 1.-PG_ge
        #
        Err_g = +1-F_e_P-PG_ge
        Err_e = 1-F_g_P -PG_eg
        #print 'errors without overlap:', Err_e, Err_g #!V
        result_dict = {'F':F_RO, 'F_g':F_g, 'F_e':F_e, 'F_post':F_RO_post,
        'F_post_g':F_g_P, 'F_post_e':F_e_P, 'F_gaus':1.-PG_eg-PG_ge, 'F_gaus_eg':1.-PG_eg,
        'F_gaus_ge':1.-PG_ge, 'Err_e':Err_e, 'Err_g':Err_g}
        return result_dict

def get_fidelity_nshot(parameters, n=3, Tabor_loading = 1, counts = 10000, average = 12):
    F_RO = 0
    F_g  = 0
    F_e  = 0
    F_RO_post = 0
    F_g_P =0
    F_e_P =0
    F_gaus   = 0
    F_gaus_eg = 0
    F_gaus_ge = 0
    Err_e = 0
    Err_g = 0

    f = get_fidelity(parameters,  Tabor_loading = Tabor_loading, counts = counts, average = average)
    #print f

    F_RO = f['F']
    F_g  = f['F_g']
    F_e  = f['F_e']
    F_RO_post = f['F_post']
    F_g_P =f['F_post_g']
    F_e_P =f['F_post_e']
    F_gaus   = f['F_gaus']
    F_gaus_eg = f['F_gaus_eg']
    F_gaus_ge = f['F_gaus_ge']
    Err_e = f['Err_e']
    Err_g = f['Err_g']

    if n == 0:
        print 'error:nShots n should not be zero'
        result_dict = {'F':F_RO, 'F_g':F_g, 'F_e':F_e, 'F_post':F_RO_post,
        'F_post_g':F_g_P, 'F_post_e':F_e_P, 'F_gaus':F_gaus, 'F_gaus_eg':F_gaus_eg,
        'F_gaus_ge':F_gaus_ge, 'Err_e':Err_e, 'Err_g':Err_g}
        return result_dict

    for i in np.arange(1, n, 1):
        f = get_fidelity(parameters,  Tabor_loading = 0, counts = counts, average = average)
        #print f

        F_RO += f['F']
        F_g  += f['F_g']
        F_e  += f['F_e']
        F_RO_post += f['F_post']
        F_g_P +=f['F_post_g']
        F_e_P +=f['F_post_e']
        F_gaus   += f['F_gaus']
        F_gaus_eg += f['F_gaus_eg']
        F_gaus_ge += f['F_gaus_ge']
        Err_e += f['Err_e']
        Err_g += f['Err_g']

    F_RO = F_RO/n
    F_g  = F_g/n
    F_e  = F_e/n
    F_RO_post = F_RO_post/n
    F_g_P =F_g_P/n
    F_e_P =F_e_P/n
    F_gaus   = F_gaus/n
    F_gaus_eg = F_gaus_eg/n
    F_gaus_ge = F_gaus_ge/n
    Err_e = Err_e/n
    Err_g = Err_g/n

    result_dict = {'F':F_RO, 'F_g':F_g, 'F_e':F_e, 'F_post':F_RO_post,
    'F_post_g':F_g_P, 'F_post_e':F_e_P, 'F_gaus':F_gaus, 'F_gaus_eg':F_gaus_eg,
    'F_gaus_ge':F_gaus_ge, 'Err_e':Err_e, 'Err_g':Err_g}
    return result_dict

###########################################################
#
#
#               SEARCHING FUNCTIONS
#
#
###########################################################
#check that it works!
def find_pwr2(parameters, tpi_wanted, precision = 0.1, max_iteration = 6, averaging=3e3, first_tl=1):
    '''
    loking power2 ==x corresponded to T-pi = tpi_wanted ==y
    using Secant method. returns power2.
    IMPROVE: stop criteria, different averaging, errors
    tpi_wanted in seconds, inside transforms to nanoseconds and back
    '''
    #precision given in nanoseconds - convering to nanoseconds
    precision_sec = precision*1e-9

    if precision<0:
        print 'precision should be >0'
        return None
    if (tpi_wanted < 10*1e-9) or (tpi_wanted > 100*1e-9):
        print 'tpi is not in range 10-100ns'
        return None
    if parameters.power2 > 0:
        parameters.power2 = -0.0
        print 'power2 had wrong sign. Changed to -0.0'

    #y1,y2,y3 in nanoseconds
    x1 = parameters.power2
    if COMMENTS:
        print 'power2:',x1,'start Rabi...'
    y1 = get_rabi_pi(parameters, x1, averaging=2e3, Tabor_loading = first_tl)
    if COMMENTS:
        print 'x1=',x1,'y1=',y1

    x2 = parameters.power2*1.1 - 0.3
    if COMMENTS:
        print 'power2:',x2,'start Rabi...'
    y2 = get_rabi_pi(parameters, x2, averaging=2e3, Tabor_loading=0)
    if COMMENTS:
        print 'x2=',x2,'y2=',y2

    for i in range(max_iteration):
        x3 = (x2-x1)/(y2-y1)*(tpi_wanted-y1)+x1
        last = False
        if x3 < -10:
            print 'WARNING: rabi power2 more than -10. Bad regime'
            last = True
            x3 = -10.0
        elif x3 >0:
            print 'WARNING: rabi power2 less than 0. Bad regime'
            last = True
            x3 = -0.0
        if COMMENTS:
            print 'power2:', x3

        y3 = get_rabi_pi(parameters, x3, Tabor_loading=0, averaging=averaging)
        if COMMENTS:
            print  'Iteration:', i, '\nx3=', x3, 'tpi', y3
        if (math.fabs(tpi_wanted-y3) < precision_sec) and (math.fabs(tpi_wanted-y2)):
            return x3
        x1 = x2
        x2 = x3
        y1 = y2
        y2 = y3
        if last == True:
            print 'value on the edge returned'
            return x3

    print  'Warning! [power2,Tpi] is not precise enough! \n', 'power2=',x3,'tpi=',y3
    return x3

def find_pwr2_old(parameters, tpi_wanted, precision = 0.1, max_iteration = 6):
    '''
    loking power2 ==x corresponded to T-pi = tpi_wanted ==y
    using Secant method. returns power2.
    IMPROVE: stop criteria, different averaging, errors
    tpi_wanted in seconds, inside transforms to nanoseconds and back
    '''
    #transform from seconds to nanoseconds
    tpi_wanted = 1e9*tpi_wanted

    if precision<0:
        print 'precision should be >0'
        return None
    if (tpi_wanted < 10) or (tpi_wanted > 100):
        print 'tpi is not in range 10-100ns'
        return None
    if parameters.power2 > 0:
        parameters.power2 = -0.0
        print 'power2 had wrong sign. Changed to -0.0'

    #y1,y2,y3 in nanoseconds
    x1 = parameters.power2
    if COMMENTS:
        print 'power2:',x1,'start Rabi...'
    y1 = 1e9*get_rabi_pi(parameters, x1, averaging=2e3)
    if COMMENTS:
        print 'x1=',x1,'y1=',y1

    x2 = parameters.power2*1.1 - 0.3
    if COMMENTS:
        print 'power2:',x2,'start Rabi...'
    y2 = 1e9*get_rabi_pi(parameters, x2, averaging=2e3, Tabor_loading=0)
    if COMMENTS:
        print 'x2=',x2,'y2=',y2

    for i in range(max_iteration):
        x3 = (x2-x1)/(y2-y1)*(tpi_wanted-y1)+x1

        if x3 < -10:
            print 'power2 more than -10. Bad regime'
            return None
        elif x3 >0:
            print 'power2 less than 0. Bad regime'
            x3 = -0.0
            print 'power2=-0.0'
            return None
        if COMMENTS:
            print 'power2:', x3
        y3 = 1e9*get_rabi_pi(parameters, x3, Tabor_loading=0, averaging=3e3)
        if COMMENTS:
            print  'Iteration:', i, '\nx3=', x3, 'tpi', y3
        if (math.fabs(tpi_wanted-y3) < precision) and (math.fabs(tpi_wanted-y2)):
            return x3
        x1 = x2
        x2 = x3
        y1 = y2
        y2 = y3

    print  'Warning! [power2,Tpi] is not precise enough! \n', 'power2=',x2,'tpi=',y2
    return x3

def find_fqubit(parameters, precision = 1, max_iteration = 5, averaging = 1e3):
    df = 0

    for i in range(max_iteration):
        df = get_ramsey_df(parameters, averaging = averaging)
        parameters.freq_q = parameters.freq_q - df*1e-3
        if df < precision:
            return parameters.freq_q
        aver = aver*2

    print 'not success: df= ', df, 'MHz'
    return None

def get_closest_tpi_old(params):
    '''
    Function check tpi on given pwr2 and find closest even tpi
    With returned tpi one could run find_pwr2(param, tpi) to find pwr2
    '''
    tpi_raw = get_rabi_pi(param0, param0.power2)
    if COMMENTS:
        print 'tpi_raw',tpi_raw
    tpi_wanted = 1e-9*round(tpi_raw*1e9)
    if COMMENTS:
        print 'after round', tpi_wanted
    if (round(tpi_wanted*1e9) % 2) != 0:
        tpi_wanted = tpi_wanted + 1e-9
    if COMMENTS:
        print 'after check even:', tpi_wanted
    return tpi_wanted

def get_closest_tpi(params):
    '''
    Function check tpi on given pwr2 and find closest even tpi
    With returned tpi one could run find_pwr2(param, tpi) to find pwr2
    Use function 'closest_even()'
    Improve: if closest even too small for achive?
    '''
    tpi_raw = get_rabi_pi(param0, param0.power2)
    if COMMENTS:
        print 'tpi_raw',tpi_raw
    tpi_wanted = 1e-9*closest_even(tpi_raw*1e9)
    if COMMENTS:
        print 'after closest_even()', tpi_wanted
    return tpi_wanted

def set_closest_tpi(params):
    '''
    Do the Rabi, find closest even tpi, find power2 for this tpi
    Set this values to params
    return [power2, tpi]
    '''
    if COMMENTS:
        print 'start to get closest tpi'

    tpi_wanted = get_closest_tpi(params)

    if COMMENTS:
        print 'closest tpi got', tpi_wanted

    pwr2 = find_pwr2(params, tpi_wanted, precision = 0.05)

    if COMMENTS:
        print 'pwr2 got'

    params.tpi = tpi_wanted
    params.power2 = pwr2
    return [pwr2, tpi_wanted]
