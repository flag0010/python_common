from collections import deque
mean = lambda s: sum(s)*len(s)**-1

def sliding_window(x, window_size, output_interval=10, flip_start=False, summary_stat=mean, start_output_frac = 1.0):
    '''simple sliding window implementation
       x = input values you want to slide windows over
       window_size = size of window
       output_interval = the number of strides you want the window to take before outputing a value, min=1, max=window_size
       flip_start = logical, do you wan to flip the first window and tack it onto the beginning.
                    i.e. if x = [0,1,2,3,4,5,6,7,8] and window size is 3, then with flip_start=True
                    x becomes [2,1,0,0,1,2,3,4,5,6,7,8]. This fills the sliding window with a mirror image of the beginning values
                    and is a decent way (sometimes) to avoid outputing nothing while the window is filling
       summary_stat = any func that calcs a summary stat on the window, default=mean, but you could use median, sd, etc.
       start_outputing_frac = how many values in the window before you start outputing values, default=1.0 meaning wait
                    till window_size values are in before outputting anything.  Must be between 0-1'''
    start_output_frac = float(start_output_frac)
    if flip_start:
        x = list(reversed(x[:window_size])) + x
    d = deque(maxlen=window_size)
    idx = 0
    for i in x:
        d.append(i)
        idx+=1
        if idx == output_interval and len(d) >= window_size * start_output_frac:
            yield summary_stat(d)
        if idx == output_interval: idx = 0


if __name__ == '__main__':
    import random
    from math import sin
    from matplotlib import pyplot as plt
    k = [random.random()+sin(i/100.) for i in xrange(10000)]
    x = list(sliding_window([i/100. for i in range(10000)], 10))
    y = [i for i in sliding_window(k, 10)]
    plt.scatter([i/100. for i in range(10000)], k, alpha=0.1, marker='.')
    plt.plot(x,y)
    plt.show()