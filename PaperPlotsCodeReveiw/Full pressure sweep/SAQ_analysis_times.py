#!/usr/bin/env python3

import os
import sys
import ROOT
import numpy as np

from SAQ_DAQ import N_SAQ_CHANNELS


def fix_wrap_around(times):
    """
    update timestamp values from the zybo
    every timestamp that is less than the previous indicates a 'wrap around'
    """

    wrap_val = 2**32
    n_wrap = 0
    inc_times = [times[0]]
    for i, t in enumerate(times[:-1]):

        if times[i+1] < t:
            n_wrap += 1

        cur_time = times[i+1] + (n_wrap * wrap_val)
        inc_times.append(cur_time)

    return inc_times

def filter_saq(resets, SAQ_DIV, ZYBO_FRQ, setXtoTime, noiseFil, min_time=10e-6*2, min_binBit=0):
    """
    Ensure that resets time differences exclude impossible minimum time from
    length of reset pulse.
    ARGS:
       resets       : list of 32bit timestamp values from Zybo
       SAQ_DIV      : clock division register on the zybo. This should be read in from metadata
       ZYBO_FRQ     : zybo norminal frequecy in hertz. This should be read in from metadata
       min_time     : time in seconds of reset pulse width, default is 10us from SAQ
       min_binBit   :
       noiseFil     : boolean, determines whether noise is filtered
    RETURNS:
       rtds         : parsed list of reset time differences
    """
    assert isinstance(resets, list), "expect a list of resets to act on"
    # build the reset list with RTDs
    rtds = []
    iLastReset = 0
    for i, r in enumerate(resets[:-1]):
        rtd = resets[i + 1] - resets[iLastReset]
        if rtd < 0:
            rtd += 2 ** 32
        if noiseFil ==True:
            if (rtd * SAQ_DIV) / ZYBO_FRQ > min_time and rtd > min_binBit:
                if setXtoTime == True:
                    rtds.append((rtd * SAQ_DIV) / ZYBO_FRQ)
                if setXtoTime == False:
                    rtds.append(rtd)
               
                iLastReset = i + 1
        else:
            if (rtd * SAQ_DIV) / ZYBO_FRQ > min_time:
                print(min_time)
                if setXtoTime == True:
                    #print("tru")
                    #print("timeConv: " + str((rtd * SAQ_DIV) / ZYBO_FRQ) + "rtdOriginal: " +str(rtd) + "saq_div: " + str(SAQ_DIV)+ "ZYBO_FRQ: " + str(ZYBO_FRQ)  )
                    rtds.append((rtd * SAQ_DIV)/ ZYBO_FRQ)
                if setXtoTime == False:
                    rtds.append(rtd)
                    print("False")

               	iLastReset = i + 1

    iLastReset = 0
    if not resets:
    	times = [0]
    else:
        times = [resets[iLastReset]]
    for i, r in enumerate(resets[:-1]):
    	times.append(resets[i+1])
    
    return rtds, fix_wrap_around(times)


def hist_rtd(rtd, gaussCh, setXtoTime, ch, i=0):
    """
    Helper function which can take in a list of an RTD (usually output of
    filter_saq method) and return a filled TH1F with an applied gaus fit.

    This histogram method is intended to find the mean of the fit for drift
    current, therefore the algorithm we use will be to find the mean of the list
    of rtd, the std of the values, and we will supply bins of of sqrt(len(rtd))
    with a bin width of 3 sigma before applying a fit. This helps the fitter
    find the region of interest (ROI).
    ARGS:
       rtd      : list of RTD in purely timestamp format, conversion into *time*
             should be *last* step
       ch       : number of channel for histogram, which stores histogram name
       gausCh   : boolean, True for gaus fit, False for standard
    RETURNS:
       h        : TH1F histogram
       mu       : mean value of gaus fit
       sig      : sigma value of returned gaus fit
    """
    print("rtd: " + str(len(rtd)))
    if gaussCh== True:
        # if no resets, just make a histogram
        # or len(rtd) == 1
        
        if len(rtd) == 0:
            nbins = int(np.sqrt(len(rtd))) + 100
            h = ROOT.TH1F(f"h_{ch}", "Ch:" + str(ch + 1) + " (Avg Current(A): 0)", nbins, 0, 0)
            for r in rtd:
                h.Fill(r)
            xaxis = h.GetXaxis()
            xaxis.SetTitle("Delta T (s)")
            xaxis.SetLabelSize(0.03)

            yaxis = h.GetYaxis()
            yaxis.SetTitle("# of Events")
            yaxis.SetLabelSize(0.03)

            h.Draw()
            return h

        mu, sig = 0, 0

        # build parameters of TH1F
        mean, std = np.mean(rtd), np.std(rtd)
        #avgCurrent = 10e-12 * 4 * 0.015 / mean
        avgCurrent = 10e-12 * 3.4 / mean
        nbins = int(np.sqrt(len(rtd))) + 100
        h = ROOT.TH1F(f"h_{ch}", "Ch:" + str(ch + 1) + " (Avg Current(A):" + str("{:.1e}".format(avgCurrent)) + ")",
                      nbins, 0, mean + 3 * std)
        for r in rtd:
            h.Fill(r)
        xaxis = h.GetXaxis()
        if setXtoTime == False:
            xaxis.SetTitle("Delta T (cycles)")
        if setXtoTime == True:
            xaxis.SetTitle("Delta T (s)")
        xaxis.SetLabelSize(0.03)

        yaxis = h.GetYaxis()
        yaxis.SetTitle("# of Events")
        yaxis.SetLabelSize(0.03)

        gaus = ROOT.TF1("gaus", "gaus", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        gaus.SetParameters(h.GetMaximum(), h.GetMean(), h.GetRMS())
        gaus.SetParLimits(1, h.GetMean() - 3 * h.GetRMS(), h.GetMean() + 3 * h.GetRMS())
        gaus.SetParLimits(2, 0, 2 * h.GetRMS())

        h.Fit("gaus", "RQ")
        h.Draw()
        gaus.Draw("same")
        
        
        if len(rtd) < 3:
            sig = 0
            mu=0
        else:
            sig = h.GetFunction("gaus").GetParameter("Sigma")
            mu = h.GetFunction("gaus").GetParameter("Mean")
      
        
        return h, mu, sig
    else:
        # if no resets, just make a histogram
        if len(rtd) == 0:
            nbins = int(np.sqrt(len(rtd))) + 100
            h = ROOT.TH1F(f"h_{ch}", "Ch:" + str(ch + 1) + " (Avg Current(A): 0)", nbins, 0, 0)
            for r in rtd:
                h.Fill(r)
            xaxis = h.GetXaxis()
            xaxis.SetTitle("Delta T (s)")
            xaxis.SetLabelSize(0.03)

            yaxis = h.GetYaxis()
            yaxis.SetTitle("# of Events")
            yaxis.SetLabelSize(0.03)

            h.Draw()
            return h

        mu, sig = 0, 0

        # build parameters of TH1F
        mean, std = np.mean(rtd), np.std(rtd)
        #avgCurrent = 10e-12 * 4 * 0.015 / mean
        avgCurrent = 10e-12 * 3.4 / mean
        nbins = int(np.sqrt(len(rtd))) + 100
        h = ROOT.TH1F(f"h_{ch}", "Ch:" + str(ch + 1) + " (Avg Current(A):" + str("{:.1e}".format(avgCurrent)) + ")",
                      nbins,
                      0, mean + 3 * std)
        for r in rtd:
            h.Fill(r)
        xaxis = h.GetXaxis()
        xaxis.SetTitle("Delta T (s)")
        xaxis.SetLabelSize(0.03)

        yaxis = h.GetYaxis()
        yaxis.SetTitle("# of Events")
        yaxis.SetLabelSize(0.03)

        h.Draw()

        return h, mean, std


def main(input_file, use_multithread=True):
    """
    Test script for running a simple analysis on ROOT files generated from
    qdb_interface based GUIs like SAQ_DAQ.py

    Supply input file to read, then this script should be able to parse that
    data and supply desired graphics
    """
    if use_multithread:
        ROOT.EnableImplicitMT()  # gotta go fast

    # open up ttree into an rdataframe, which is easy to convert to a list
    rdf = ROOT.RDataFrame("tt", input_file)

    # rip everything immediately into a dictionary, where the keys are
    # the branch names
    data = rdf.AsNumpy()

    # numpy arrays of the data we need
    ts = data["Timestamp"]
    masks = data["ChMask"]

    # snag the meta data from the tfile
    # these values are defined in make_root.py script when ROOT file is created
    # from binary
    tf = ROOT.TFile(input_file, "READ")
    meta_data = tf.mt
    SAQ_DIV = -1
    ZYBO_FRQ = -1
    version = 0

    # Grab these settings from the tf file
    for evt in meta_data:
        SAQ_DIV = evt.SAQ_DIV
        version = evt.Version
        ZYBO_FRQ = evt.Zybo_FRQ
    assert version >= 0x3f, f"version of root file is too old! 0x{version:02x} <= 0x3f"
    # assert SAQ_DIV >= 1, f"SAQ_DIV not properly defined: {SAQ_DIV} not >= 1"
   
    assert ZYBO_FRQ >= 30e6, f"ZYBO_FRQ not properly defined: {ZYBO_FRQ} not >= 30 MHz"

    # make a quick way to ensure the channel we want is in the mask
    m = lambda ch, mask: 1 << ch & mask

    # create a list of the channels and all of their resets
    chResets = [[t for t, mask in zip(ts, masks) if m(ch, mask)] for ch in range(N_SAQ_CHANNELS)]

    # inspect the resets to ensure they make sense
    for ch, resets in enumerate(chResets):
        print(f"ch: {ch + 1} has {len(resets)} resets.")

    # chRTD is a list of a list where the first index is the channel number -1
    # (lists are zero counted), which contain the sequential time since last
    # reset data

    ###################################
    # extra analysis can proceed here #
    ###################################
    setXtoTime = True
    noiseFil = False
    chRTDTimes = [filter_saq(r, SAQ_DIV, ZYBO_FRQ, setXtoTime, noiseFil) for r in chResets]

    chRTD = [rtds for rtds, times in chRTDTimes]
    chTimes = [times for rtds, times in chRTDTimes]

    # let's make some gaussian histograms of all of the RTDs we have
    outf = ROOT.TFile("saqAna.root", "RECREATE")
    hists = []
    gaussCh = True
    
    print("chrtd: "+str(len(chRTD)))

    # Open a .txt file to save the rtds to for processing
    file = open(input_file.split(".root")[0]+'.txt', 'w')
    file.write(str(chRTD))
    file.close()

    file = open(input_file.split(".root")[0]+'_times.txt', 'w')
    file.write(str(chTimes))
    file.close()

    # Make a 16 panel root histogram with the data
    for ch, rtd in enumerate(chRTD):

        # if no resets, make graphs but leave them empty
        if len(rtd) == 0:
            h = hist_rtd(rtd, gaussCh, setXtoTime, ch)  
            hists.append(h)
            h.Write()
            continue

        # Otherwise make a nice plot
        h, mean, sig= hist_rtd(rtd, gaussCh, setXtoTime, ch)
        hists.append(h)

        # time conversion lambda
        t = lambda x: x * SAQ_DIV / ZYBO_FRQ

        # print and save results
        print(f"Channel-{ch} has mean={t(mean):2e} and sigma={t(sig):2e}")
        h.Write()

    canvas = ROOT.TCanvas("canvas", "Histograms", 3500, 3500)
    canvas.SetCanvasSize(3500, 3500)
    canvas.Divide(4, 4)

    for i, hist in enumerate(hists):
        canvas.cd(i + 1)
        hist.Draw()
        hist.Write()
        canvas.Modified()
        canvas.Update()
    canvas.Write()
    
    canvas.SaveAs(input_file.split(".root")[0]+".png")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("ERROR user must supply input file name")
    else:
        input_file = sys.argv[1]
        if not os.path.isfile(input_file) or os.path.getsize(input_file) == 0:
            print("empty or non-existent input ROOT file", input_file)
            sys.exit(-1)
        main(input_file)
