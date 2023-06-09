def uta_area():
    #Contains the start,end of each copper ring consecuatively
    start_and_end_Channel = [0,0.53,0.79,1.3,1.47,2.01,2.19,2.89,3.08,3.93,4.1,4.9,5.1,5.9,6.12,7.37,7.58,9.92,10.11,12.39,12.59,14.87,15.090,19.8,20.12,24.79,25.12,29.86,30.15,39.76,40.060,47.7]

    #split into the start and end of each channel
    end_ch, start_ch = start_and_end_Channel[1::2], start_and_end_Channel[::2]

    #make a 2d list of min and max for each channel, where
    # the start and end is the middle of the space between each copper ring
    minmax=[]

    #start at 0mm
    end=0
    num=0
    for i in range (0,16):
        num+=1
        start=end
        if num==16:
            #If the last channel, count the end as the end
            end = round(end_ch[i], 3)
        else:
            #If not the last channel,count the end as the (end+half the distance to the next ring)
            end=round(end_ch[i]+((start_ch[i+1]-end_ch[i])/2),3) #round to three decimal places
        startend=[start,end]
        minmax.append(startend)
        ring_areas = [round((3.14159 * (end ** 2)) - (3.14159 * (start ** 2)), 3) for start, end in annuli]
    return ring_areas

def weseley_area():
    annuli = [[ 0.000,  0.670], [ 0.670,  1.386], [ 1.386,  2.106], [ 2.106,  2.981],
           [ 2.981,  4.005], [ 4.005,  4.994], [ 4.994,  6.015], [ 6.015,  7.481],
           [ 7.481,  9.994], [ 9.994, 12.497], [12.497, 15.023], [15.023, 19.996],
           [19.996, 24.962], [24.962, 30.026], [30.026, 39.977], [39.977, 50.065]]

    ring_areas = [round((3.14159 * (end**2))-(3.14159*(start**2)), 3) for start, end in annuli]
    return ring_areas

print(weseley_area())
