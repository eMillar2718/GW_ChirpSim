import bilby
import numpy as np
from gwpy.timeseries import TimeSeries
import os
import random
import subprocess

from gravityspy.plot.plot import plot_qtransform



def get_gwosc_data(number_of_injections, sample_rate):

    #Specifying GPS time ranges for each observing run
    O1_range = [1126051217, 1137254417]
    O2_range = [1164556817, 1187733618]
    O3a_range = [1238166018, 1253977218]
    O3b_range = [1256655618, 1269363618]

    Obs_runs = [O1_range, O2_range, O3a_range, O3b_range]

    merger_times = []
    start_times = []
    end_times = []
    time_series_all = []

    for i in range(int(number_of_injections)):

        for attempt in range(10): #Try generated GPS time 10 times

            #Randomly select an observing run
            run_select = Obs_runs[random.randint(0,len(Obs_runs)-1)]

            print('Run Selection = {}'.format(run_select))

            #Randomly select a start time within that observing run, leaving room at the end for signal trail-off
            #Specify start times, end times, merger times

            merger = random.randint(run_select[0], run_select[1]-duration)
            start = merger - duration/2
            end = start + duration

            try: #Attempt to locate GWOSC dataset 
                time_series = TimeSeries.fetch_open_data('H1', start, end, sample_rate = sample_rate, verbose = True)
            
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
            time_series = TimeSeries.fetch_open_data('H1', start, end, sample_rate = sample_rate, verbose = True) #known working value
            
        time_series_all.append(time_series)

    return time_series_all, merger_times, start_times, end_times

def simulate_and_inject_waveforms(merger_times, timeseries_dataset): 

    """

    merger_times: array of gps times to specify as the merger time

    timeseries_dataset: gwosc timeseries data to inject signals into
    
    """



    injection_parameters_all = []
    waveform_arguments_all = []
    waveform_generators_all = []

    for i in range(len(timeseries_dataset)):

        mass_1_gen = np.random.randint(5,50)
        mass_2_gen = np.random.randint(5,mass_1_gen)

        injection_parameters_single = dict(
            mass_1=mass_1_gen, 
            mass_2=mass_2_gen,
            a_1=0.4, # spin
            a_2=0.3, # spin
            tilt_1=0.5, # spin
            tilt_2=1.0, # spin
            phi_12=1.7, # spin
            phi_jl=0.3, # spin
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

    for i in range(len(timeseries_dataset)):

        data_injected, metadata = bilby.gw.detector.inject_signal_into_gwpy_timeseries(data = timeseries_dataset[i],
            waveform_generator=waveform_generators_all[i], parameters=injection_parameters_all[i], det = 'H1'
        )

        injections.append(data_injected)
        metadatas.append(metadata)

    return injections, metadatas

def save_to_gwf(data, start_times, file_path):

    if not os.path.exists(file_path):
        os.mkdir(file_path)
        
    for i in range(len(data)):

        # Naming gwf channel to the Omicron specification
        data[i].name = "H1:SIM-CHIRP"
        data[i].channel = "H1:SIM-CHIRP"

        # Writing injection as .gwf file
        data[i].write(file_path + '/H-H1_SIM-{0}-{1}.gwf'.format(int(start_times[i]), duration))

def pass_injection_through_omicron(start_times, end_times, duration, file_path = './injections'):

    """
    injections: list of timeseries data with injected signal
    start_times: list of of the GPS start time for each timeseries data
    duration: duration of the timeseries
    
    """
    snr_all = []
    frequency_all = []
    peak_time_all = []

    for i in range(int(len(start_times))):

        #Generating Cache file, points to the location of the injection

        omicron_cache_file = open('./sim.lcf', "w+")
        omicron_cache_file.write('file://localhost' + cwd + file_path + '/H-H1_SIM-{0}-{1}.gwf'.format(int(start_times[i]), duration))
        omicron_cache_file.close()

        # Run injection through Omicron

        # Setting the terminal command to run omicron
        omicron_command = "omicron-process --gps {0} {1} --ifo H1 --config-file ./sim.ini --output-dir ./run --no-segdb --cache-file ./sim.lcf -vvv --file-tag SIM SIM --no-submit".format(int(start_times[i]),int(end_times[i]))

        # running pyomicron through terminal and printing the output
        run_pyomicron = subprocess.run(omicron_command.split(), shell=False, capture_output=True, text=True)
        #print(run_pyomicron.stdout)

        # running omicron bash script generated by pyomicron 
        run_omicron_script = subprocess.run("./run/condor/omicron-SIM.sh", shell=True, capture_output=True, text=True)
        #print(run_omicron_script.stdout)

        # specifying the path of the omicron trigger files
        omicron_output_path = cwd + "/run/merge/H1:SIM-CHIRP/H1-SIM_CHIRP_OMICRON-{}-124.root".format(int(start_times[i]+2)) #+2 as the chunk starts 2s later

        # printing the result of the trigger files
        print_omicron_all = subprocess.run(["omicron-print", "file={}".format(omicron_output_path)], shell=False, capture_output=True, text=True)
        print_omicron_snr = subprocess.run(["omicron-print", "file={}".format(omicron_output_path), "print-freq=0"], shell=False, capture_output=True, text=True)

        #print(print_omicron_all.stdout)
        print(print_omicron_snr.stdout)

        #Formatting the output readings

        output_split = print_omicron_all.stdout.split("\n")
        output_numbers = output_split[4:-1]

        #print(output_numbers)

        output_array = []
        for row in output_numbers:
            output_array.append(list(map(float, row.split())))

        #print(output_array)

        peak_times = []
        frequencies = []
        snrs = []

        for i, readings in enumerate(output_array):
            peak_time = output_array[i][0]
            frequency = output_array[i][1]
            snr = output_array[i][2]

            peak_times.append(peak_time)
            snrs.append(snr)
            frequencies.append(frequency)

        snr_all.append(snrs)
        frequency_all.append(frequencies)
        peak_time_all.append(peak_times)

    return peak_time_all, frequency_all, snr_all

def check_omicron_threshold(snrs, start_times, injection_path):

    """
    Moves successful injections into a /successful_injections folder and deletes the failed ones (minimum SNR of 7.5)

    snrs: list of snr values from omicron
    start_times: omicron start times
    injection_path: path to injections 

    """

    for i in range(int(len(start_times))):

        if snrs:

            print('max SNR for injection {0} is {1}'.format(i+1,max(snrs[i])))
        
        else:

            print('No trigger file detected for injenction')

    for i in range(len(start_times)):


        if max(snrs[i]) > 7.5:
            print('SNR greater than 7.5, injection allowed')

            os.mkdir('./successful_injections')
            os.rename(injection_path + '/H-H1_SIM-{0}-{1}.gwf'.format(int(start_times[i]), duration), './successful_injections' + '/H-H1_SIM-{0}-{1}.gwf'.format(int(start_times[i]), duration))

        else:
            print('SNR less than 7.5, injection rejected')
            os.remove(injection_path + '/H-H1_SIM-{0}-{1}.gwf'.format(int(start_times[i]), duration))

def generate_qscans(injections, start_times, merger_times, default_path = './plots'):    


    print('Generating Q-Scans...')

    durations = [0.5, 1, 2, 4]
    plot_normalized_energy_range = (0,25)
    plot_time_ranges = durations
    detector_name = 'H1'

    qspecgrams = {}
    ind_figs_all = {}
    superfigs = []

    for i in range(len(injections)):
        
        qspecgram_durations = []
        for j in range(len(durations)):

            print('Plot {} Generated'.format(j+1))
            
            qspecgram = injections[i].q_transform(frange = [10,2048], outseg=(merger_times[i] - durations[j]/2, merger_times[i] + durations[j]/2),tres=0.002, fres = 0.5, whiten=True, qrange = [4,64], gps = merger_times[i])

            qspecgram_durations.append(qspecgram)
            
        qspecgrams["injection {0}".format(i+1)] = qspecgram_durations

        ind_fig_all, superfig = plot_qtransform(qspecgrams["injection {}".format(i+1)], plot_normalized_energy_range, plot_time_ranges, detector_name, start_times[i])

        ind_figs_all["injection {0}".format(i+1)] = ind_fig_all
        superfigs.append(superfig)

    answer = "x"
    while answer not in ["Yes", "yes", "Y", "y", "No", "no", "N", "n"]:
        print('Use default output path for plots? "{0}"'.format(default_path))
        answer = input("[y/n] ")

        if answer == "Yes" or answer == "yes" or answer == "Y" or answer == "y":

            if not os.path.exists(default_path):
                os.mkdir(default_path)

        elif answer == "No" or answer == "no" or answer == "N" or answer == "n":

            print('Please specify output path')
            custom_path = input()

            if not os.path.exists(custom_path):
                os.mkdir(custom_path)

    for i in range(len(injections)):

        for j in range(len(durations)):

            ind_figs_all["injection {0}".format(i+1)][j].savefig(default_path + '/qscan{0}-{1}-{2}s.jpg'.format(i+1, merger_times[i], durations[j]))

    print('Q-scans generated successfully, saved to {0}'.format(default_path))

#Universal parameters
sampling_frequency = 16384
duration = 128
minimum_frequency = 10

cwd = os.getcwd()

number_input = input('Enter number of injections to generate: ')

timeseries_all, merger_times, start_times, end_times = get_gwosc_data(number_of_injections=number_input, sample_rate=sampling_frequency)

save_to_gwf(data=timeseries_all, start_times=start_times, file_path= './ligo_data')

peak_times, frequencies, snrs = pass_injection_through_omicron(start_times=start_times, end_times=end_times, duration=duration,file_path='./ligo_data')

injections, metadatas = simulate_and_inject_waveforms(merger_times=merger_times, timeseries_dataset=timeseries_all)

save_to_gwf(data=injections, start_times=start_times, file_path= './injections')

number_of_injections = int(len(injections))

print('Injection generated successfully, running through Omicron...')

peak_times, frequencies, snrs = pass_injection_through_omicron(start_times=start_times, end_times=end_times, duration=duration)

check_omicron_threshold(snrs=snrs, start_times=start_times, injection_path='./injections')