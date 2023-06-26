import subprocess
import os
import glob
import sys

def total_readout_time_dcm(dcm_file):
    """ Total readout time calculation from DICOM header

    This is a tool designed to extract the total readout time from a dicom file for further processing, e.g. with FSL.
    It relies on the function dcmdump and is designed to work with Philips dicom standards.
    
    Parameters
    ----------
    dcm_file: file
        Input DICOM-formated file (e.g. one file from the DICOM stack of interest)

    Returns
    ----------
    float:
        Calculated readout time

    Notes
    ----------
    Original script by Luka Liebrand (19-10-2017)
    """

    if not os.path.exists(dcm_file):
        print("Dicom file "+dcm_file+" does not exist!")
        sys.exit(1)

    ysize = "dcmdump "+dcm_file+"| grep NumberOfPhaseEncodingSteps | awk 'NR==1{print$3}' | tr -d '[]'"
    ysize=float(subprocess.check_output(ysize,shell=True,stderr=subprocess.STDOUT))
    # number of k-lines (# of PE-steps)     

    etl="dcmdump "+dcm_file+"| grep EchoTrainLength | awk 'NR==1{print$3}' | tr -d '[]' "
    etl=float(subprocess.check_output(etl,shell=True,stderr=subprocess.STDOUT))
    # echo train length
    
    field="dcmdump "+dcm_file+"| grep MagneticFieldStrength | awk '{print$3}' | tr -d '[]' "
    field=float(subprocess.check_output(field,shell=True,stderr=subprocess.STDOUT))
    # magnetic field strength
    
    fshiftpix="dcmdump "+dcm_file+"| grep \(2001,1022\) | awk '{print$3}'| awk -F '\\' '{print$1}' | tr -d '[]' "
    fshiftpix=float(subprocess.check_output(fshiftpix,shell=True,stderr=subprocess.STDOUT))
    # water-fat shift in pixels
    
    acceleration="dcmdump "+dcm_file+"| grep ParallelReductionFactorInPlane | awk '{print$3}' | tr -d '[]'  "
    acceleration=float(subprocess.check_output(acceleration,shell=True,stderr=subprocess.STDOUT))
    # parallel acceleration factor (sense) in plane

    gmr=42.576
    # gyromagnetic ratio (1H)
    fshiftppm=3.35
    # fat-shift ppm (see Haacke et al.)

    fshifthz=gmr*field*fshiftppm
    # fat-shift (Hz)
    bwhzpix=fshifthz / fshiftpix
    # bandwidth Hz per pixel
    totalbw=bwhzpix * ysize
    # total bandwidth
    
    print("total BW = "+str(totalbw))

    echospacing= 1.0/totalbw
    # echo spacing is the time between acquiring the center of consecutive k-lines, not considering acceleration
    effechospac= echospacing / acceleration
    # effective echo spacing taking into account parallel acceleration (e.g. sense)
    # but NOT multiband

    totreadout= echospacing * (ysize-1.0)
    # total readout time from center of first k-line center to center of last k-line
    # as defined by fsl

    print("The echo spacing is: "+str(echospacing)+" seconds.")
    print("The effective echo spacing is: "+str(effechospac)+" seconds.")
    print("The total readout time is: "+str(totreadout)+" seconds.")

    return totreadout



if __name__ == "__main__":
    total_readout_time_dcm(dcm_file=sys.argv[1])