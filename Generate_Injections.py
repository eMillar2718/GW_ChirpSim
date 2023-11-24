import bilby
import numpy as np
from gwpy.timeseries import TimeSeries
import os
import random

from gravityspy.plot.plot import plot_qtransform

#Universal parameters
sampling_frequency = 16384
duration = 128 #Required duration for Omicron
minimum_frequency = 10

output_path_injections = './injections' #default output path for injections
output_path_plots = './plots' #default output path for plots (q-scans)

number_injections = input('Enter number of injections to generate: ')

#Specifying GPS time ranges for each observing run
O1_range = [1126051217, 1137254417]
O2_range = [1164556817, 1187733618]
O3a_range = [1238166018, 1253977218]
O3b_range = [1256655618, 1269363618]
#O4a

Obs_runs = [O1_range, O2_range, O3a_range, O3b_range]

merger_times = []
start_times = []
end_times = []
time_series_all = []

for i in range(int(number_injections)):

    for attempt in range(10): #Try generated GPS time 10 times

        #Randomly select an observing run
        run_select = Obs_runs[random.randint(0,len(Obs_runs)-1)]

        print('Run Selection = {}'.format(run_select))

        #Randomly select a start time within that observing run, leaving room at the end for signal trail-off
        #Specify start times, end times, merger times

        merger = random.randint(run_select[0], run_select[1]-duration)
        start = merger - duration/2
        end = start + duration

        print(start, end)

        try: #Attempt to locate GWOSC dataset 
            time_series = TimeSeries.fetch_open_data('H1', start, end, sample_rate = sampling_frequency, verbose = True)
        
        except: 
            print('Could not find dataset, retrying')

        else: #if successful, break
            merger_times.append(merger)
            start_times.append(start)
            end_times.append(end)
            break

    else: #for-else: if the loop doesn't break (no successful dataset), print error message stating number of tries
       
        start = 1126259467
        start_times.append(start)
        
        end = 1126259497
        end_times.append(end)
        
        merger = start + duration/2
        merger_times.append(merger)

        print('Could not locate dataset after {0} tries, attempting test times'.format(attempt+1))
        time_series = TimeSeries.fetch_open_data('H1', start, end, sample_rate = sampling_frequency, verbose = True) #known working value
        
    time_series_all.append(time_series)

print('merger times:{}'.format(merger_times))
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

    #Printing masses for debugging 

    mass1 = injection_parameters_single["mass_1"]
    mass2 = injection_parameters_single["mass_2"]

    print('Signal {0} Masses: M1 = {1}, M2 = {2}'.format(i+1,mass1,mass2))

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

injections = []
metadatas = []

for i in range(int(number_injections)):

    data_injected, metadata = bilby.gw.detector.inject_signal_into_gwpy_timeseries(data = time_series_all[i],
        waveform_generator=waveform_generators_all[i], parameters=injection_parameters_all[i], det = 'H1'
    )

    injections.append(data_injected)
    metadatas.append(metadata)

assert len(injections) == int(number_injections) #Ensure that we have the correct number of injections

answer = "x"

#Setting path specification, chose default or custom path

while answer not in ["Yes", "yes", "Y", "y", "No", "no", "N", "n"]:
    print('Use default output path for injections? "{0}"'.format(output_path_injections))
    answer = input("[y/n]: ")

    if answer == "Yes" or answer == "yes" or answer == "Y" or answer == "y":

        if not os.path.exists(output_path_injections):
            os.mkdir(output_path_injections)

    else:

        print('Please specify output path')
        output_path = input()

        if not os.path.exists(output_path):
            os.mkdir(output_path)

for i in range(int(number_injections)):
    injections[i].write(output_path_injections + '/H-H1-SIM-{0}-{1}.gwf'.format(merger_times[i], duration))


print('Generated {0} injections successfully (saved to {1})'.format(len(injections), output_path_injections))





#Generate Q-Scans

print('Generating Q-Scans...')

durations = [0.5, 1, 2, 4]
plot_normalized_energy_range = (0,25)
plot_time_ranges = durations
detector_name = 'H1'

qspecgrams = {}
ind_figs_all = {}
superfigs = []

for i in range(int(number_injections)):
    
    qspecgram_durations = []
    for j in range(len(durations)):

        print(type(injections[i]))
        
        qspecgram = injections[i].q_transform(frange = [10,2048], outseg=(merger_times[i] - durations[j]/2, merger_times[i] + durations[j]/2),tres=0.002, fres = 0.5, whiten=True, qrange = [4,64], gps = merger_times[i])

        qspecgram_durations.append(qspecgram)
        
    qspecgrams["injection {0}".format(i+1)] = qspecgram_durations

    ind_fig_all, superfig = plot_qtransform(qspecgrams["injection {}".format(i+1)], plot_normalized_energy_range, plot_time_ranges, detector_name, start_times[i])

    ind_figs_all["injection {0}".format(i+1)] = ind_fig_all
    superfigs.append(superfig)

answer = "x"
while answer not in ["Yes", "yes", "Y", "y", "No", "no", "N", "n"]:
    print('Use default output path for plots? "{0}"'.format(output_path_plots))
    answer = input("[y/n] ")

    if answer == "Yes" or answer == "yes" or answer == "Y" or answer == "y":

        if not os.path.exists(output_path_plots):
            os.mkdir(output_path_plots)

    else:

        print('Please specify output path')
        output_path = input()

        if not os.path.exists(output_path):
            os.mkdir(output_path)

for i in range(int(number_injections)):

    for j in range(len(durations)):

        ind_figs_all["injection {0}".format(i+1)][j].savefig(output_path_plots + '/qscan{0}-{1}-{2}s.jpg'.format(i+1, merger_times[i], durations[j]))

print('Q-scans generated successfully, saved to {0}'.format(output_path_plots))