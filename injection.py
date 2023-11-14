import bilby
import numpy as np
from gwpy.timeseries import TimeSeries
import os

#Universal parameters
sampling_frequency = 16384
duration = 30
minimum_frequency = 10

print('Enter number of injections to generate')

number_injections = input()

#Generate a random merger time in O1
merger_times = []
for i in range(int(number_injections)):
    merger = 1132283095 #np.random.randint(1126051217, 1137254417-duration)
    #Using a predetermined GPS time for all injections as there are some time gaps in O1
    #find and construct an array of useable GPS ranges to generate a random time
    merger_times.append(merger)

print('merger times:{}'.format(merger_times))

#Setting the beginning and end of each strain

start_times = []
end_times = []
for i in range(int(number_injections)):
    start = merger_times[i] - (duration/2)
    start_times.append(start)

    end = start + duration
    end_times.append(end)

print('Signal start times:{}'.format(start_times))
print('Signal end times:{}'.format(end_times))

#Setting waveform parameters

injection_parameters_all = []
waveform_arguments_all = []
waveform_generators_all = []

for i in range(int(number_injections)):

    mass_1_gen = np.random.randint(5,50)
    mass_2_gen = np.random.randint(5,mass_1_gen)

    injection_parameters_single = dict(
        mass_1=mass_1_gen, 
        mass_2=mass_2_gen,
        a_1=0.4, 
        a_2=0.3, 
        tilt_1=0.5, 
        tilt_2=1.0, 
        phi_12=1.7, 
        phi_jl=0.3,
        luminosity_distance=500., 
        theta_jn=0.4, 
        psi=2.659,
        phase=1.3, 
        geocent_time=merger_times[i],
        ra=1.375, 
        dec=-1.2108
        )
    
    injection_parameters_all.append(injection_parameters_single)

for i in range(int(number_injections)):

    #Printing masses for debugging 

    injection = injection_parameters_all[i]

    mass1 = injection["mass_1"]
    mass2 = injection["mass_2"]

    print('Signal {0} Masses: M1 = {1}, M2 = {2}'.format(i+1,mass1,mass2))

for i in range(int(number_injections)):

    #Setting up waveform generator

    waveform_argument_single = dict(
        waveform_approximant="IMRPhenomPv2",
        reference_frequency=minimum_frequency,
        minimum_frequency=minimum_frequency,
        start_time=merger_times[i] - 2
    )

    waveform_arguments_all.append(waveform_argument_single)

    waveform_generator_single = bilby.gw.WaveformGenerator(
        duration = duration,
        sampling_frequency=sampling_frequency,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_argument_single,
    )

    waveform_generators_all.append(waveform_generator_single)

time_series_all = []

for i in range(int(number_injections)):

    time_series = TimeSeries.fetch_open_data('H1', start_times[i], end_times[i], sample_rate = sampling_frequency, verbose = True)

    time_series_all.append(time_series)

injections = []
metadata = []

for i in range(int(number_injections)):

    data_injected, metadata = bilby.gw.detector.inject_signal_into_gwpy_timeseries(data = time_series_all[i],
        waveform_generator=waveform_generators_all[i], parameters=injection_parameters_all[i], det = 'H1'
    )

    injections.append(data_injected)

assert len(injections) == int(number_injections)

path = './injections'
os.mkdir(path)

for i in range(int(number_injections)):
    injections[i].write('injections/injection{0}-{1}.txt'.format(i+1, merger_times[i]))

print('Generated {0} injections successfully'.format(len(injections)))



##improvements

#Randomise GPS times by identifying data gaps 
#Allow user to specify detector, sampling rate, etc from the inputs (perhaps a pop up window with a GUI?)