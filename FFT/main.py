import numpy as np
from scipy.io import wavfile
import matplotlib.pyplot as plt


def sample_array(Samp, n):
        """
        The quantity of the sound signal samples captured may be excessive,
        The function is introduced to reduce sample density. 
        This will facilitate an easier visualization.
        
        n : Sampling period
        """
        sampled_values = []
        for i in range(0, len(Samp), n):
            sampled_values.append(Samp[i])
        return sampled_values

def zero_pad_array(array, target_length):
    """
    Extend N to an integer multiple of 2.
    Haven't been used ever since it was born.
    Someone has had its mission done for him.
    """
    current_length = len(array)

    if current_length >= target_length:
        return array

    num_zeros = target_length - current_length
    zero_padded_array = array + [0.0] * num_zeros

    return zero_padded_array


class HandmadeFFT:
    """
    Handmade FFT
    signal : A discretized signal in the form of list or np.ndarray.
    time : The time of sampling.
    """
    def __init__(self, signal, time):
        self.signal = signal
        self.N = len(signal)
        self.time = time
        self.ext_N = 2**(len(signal) - 1).bit_length()
        self.ext_signal = list(self.signal) + [0.0] * (self.ext_N - len(self.signal))


    def fft(self):
        if self.ext_N <= 1:
            return self.ext_signal

        even = HandmadeFFT(self.ext_signal[::2],self.time).fft()
        odd = HandmadeFFT(self.ext_signal[1::2],self.time).fft()

        T = [np.exp(-2j * np.pi * k / self.ext_N) * odd[k] for k in range(self.ext_N // 2)]

        return [even[k] + T[k] for k in range(self.ext_N // 2)] + \
               [even[k] - T[k] for k in range(self.ext_N // 2)]

    def plott(self):
        """
        Plot the waveform.
        """
        P = np.array(self.signal)
        t = np.linspace(0, self.time, self.N)
        
        plt.figure()
        plt.plot(t, P)
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
        plt.title('Waveform of the signal')
        plt.show()


    def plotf(self):
        """
        Plot the frequency spectrum.
        """
        fft_vals = np.array(self.fft())
        P2 = np.abs(fft_vals / self.ext_N)
        P1 = P2[:self.ext_N // 2] * 2
        f = np.linspace(0, (len(self.ext_signal)/self.time) / 2, self.ext_N // 2)
        
        plt.figure()
        plt.plot(f, P1)
        plt.xlabel('Frequency')
        plt.ylabel('Amplitude')
        plt.title('FFT of the signal')
        plt.show()


if __name__ == '__main__':
    """
    Several tests for the class.
    """
#################################
    """
    A rectangular signal for test.
    """
    Fs = 1000
    T = 1 / Fs
    duration = 5.0
    N = int(Fs * duration)
    rec_t = np.arange(0, duration, T)
    pulse_width = 0.1
    amplitude = 1.0
    
    signal = amplitude * np.heaviside(np.mod(rec_t, pulse_width * 2) - pulse_width, 1)
    result = HandmadeFFT(signal, 5.0)
    result.plotf()

#################################
    """
    'test1.wav' lasts for 5 sec and is a pure sound mixxed with random noise.
    """
    sample_rate, samples_1 = wavfile.read('test1.wav')
    # noise = 1000 * np.random.randn(len(samples_1),2)
    # samples_1 = np.array(samples_1) + noise

    cut_samples_1 = sample_array(samples_1, 5)
    result_wav1 = HandmadeFFT(cut_samples_1, 5.0)
    result_wav1.plotf()

#################################
    """
    'test2.wav' lasts for 20 sec and is a piece of record from 3b1b.
    """
    sample_rate, samples_2 = wavfile.read('test2.wav')
    cut_samples_2 = sample_array(samples_2, 5)
    result_wav2 = HandmadeFFT(cut_samples_2, 20.0)
    result_wav2.plotf()
