##############################################
##          Phasor Estimation Algorithms    ##
##          OpenPMU Project                 ##
##          Final Year Project              ##
##          Keith Houston                   ##
##          40026136                        ##
##          khouston03@qub.ac.uk            ##
###                                        ###
##############################################

##Functions:
##
##WrapDifference - Converts phase angle difference between +-pi range
##GenPredefinedCurve - Generate expected power signal for input frequency
##PayloadConvert - Decodes Base64
##Algorithm_CSV - Main Function, performs algorithm for CSV input/ouput
##Algorithm_UDP - Main Function, performs algorithm for UDP input/ouput


from numpy import *
from math import pi,atan,atan2,degrees,sqrt,fabs
import csv

import socket
import xml.etree.ElementTree as ET
import base64
from operator import add


def WrapDifference(Theta_Rad,Prev_Theta_Rad):
##        #WrappedOutput = (difference + (pi/2))%(pi) - (pi/2)


        #Diff_1
        Diff_1=Theta_Rad - Prev_Theta_Rad

        #Diff_2
        if(Theta_Rad>0):
                Diff_2_Present = (pi)-Theta_Rad         #Positive Theta to +pi
        else:
                Diff_2_Present = Theta_Rad-(-pi)        #Negative Theta to -pi
        if(Prev_Theta_Rad>0):
                Diff_2_Previous = (pi)-Prev_Theta_Rad   #Positive Prev_Theta to +pi
        else:
                Diff_2_Previous = Prev_Theta_Rad-(-pi)  #Negative Prev_Theta to -pi

        Diff_2 = Diff_2_Present+Diff_2_Previous         #Both should be positive difference (no direction indicated yet)

        SmallestTest = min(fabs(Diff_1),Diff_2)   #Find smallest difference

        
        if(SmallestTest == Diff_2):
                if(Theta_Rad>0) and (Prev_Theta_Rad<0): #Determine Direction
                        Smallest = Diff_2*-1
                else: Smallest = Diff_2
        else:
                Smallest = Diff_1 #Direction in itself

              
        return Smallest



def GenPredefinedCurve(m,n,Fs,Freq):
        # m - Syncrophasor computation window N = 2m+1 (typical = 128)
        # n - number of harmonics to include (typical = 50)
        # Fs - Sampling Frequency
        # Freq - Fundamental Frequency f - in absence of of previous frequency estimation (first window) use nominal system freq
        grid_m,grid_n =mgrid[-m:(m+1):1, 1:(n+1):1]               #generate grid of n and m values
                                                                  #not inclusive, 1:51 is len=50

        T=(1.0/Fs)

        RealComp = cos(2*pi*Freq*grid_m*grid_n*T)                 # columns of Cos(-nwmT)
        ImajComp = -1*sin(2*pi*Freq*grid_m*grid_n*T)              # columns of -Sin(-nwmT)

        idx = arange(1,(n+1),1)                                 #Indexing of merging Real and Imaj Components
        RealandIm = insert(RealComp,idx,ImajComp,axis=1)        #Join matrix together


        PredefinedCurve =zeros(((2*m+1),(2*n+2)))         #Create empty array
                                                          #Size of RealandIm + 2 columns for DC component and m


        #Input Real&Imaj components, column of 1 and column of m values into the Zero Array
        PredefinedCurve[:,:-2]=RealandIm
        PredefinedCurve[:,(2*n)]=1
        PredefinedCurve[:,(2*n+1)]=arange((-m),(m+1),1)

        return PredefinedCurve


def PayloadConvert(Payload_base64):
    ##this conversiion is for two 4-hex per number, direct b64decode may be used for 8-hex numbers
    ##[::2] and {1::2] collect every other item starting with first or second respectively.

    Payload_4hexDec = [ord(c) for c in base64.b64decode(Payload_base64)] #convert all 512 4-byte hex to decimal
    MSB = [x*256 for x in Payload_4hexDec[::2]] #multiply MSB by 256
    Payload_Output = map(add,MSB,Payload_4hexDec[1::2]) #add LSB . Output = ((MSB*256+LSB)) output is 256 integers
    return Payload_Output



def Algorithm_CSV(Freq,Fs,NumHarmonics,FileName):

    
    Data_List = genfromtxt(FileName,delimiter=',')      #Buffer input file to list
    csvwrite = csv.writer(open("Phasor_"+FileName,"wb"))#identify csv to write phasor results to.


    Samples_Buffer = 257                                        #Size of computational window
    m = (Samples_Buffer-1)/2                                    #Phasor estimated for middle of computational window
    StartOffset = m                                             #Put centre of computational window at start of 2nd cycle
    k=(Samples_Buffer-1)                                        #Amount of samples windows is moved each phasor

    CompCycles = int(len(Data_List)/Samples_Buffer)-1  #Number of runs to cover all the data.


    for x in range(0,CompCycles):             #Loop for number of phasors to calculate
            
        TimeStamp=(float(x+1)/50)       #(x+1) as first waveform is offset, estimation starts at start of waveform.
        if(x==0):
            CompWindowBuffer=Data_List[(0+StartOffset):(Samples_Buffer+StartOffset)] #take number of cycles in buffer
        else:
                CompWindowBuffer=Data_List[((Samples_Buffer-1)*x+StartOffset):((Samples_Buffer-1)*x+Samples_Buffer+StartOffset)]


        for i in range(0,2):
                #print m,NumHarmonics,Fs,Freq
                A = GenPredefinedCurve(m,NumHarmonics,Fs,Freq)       #Generate predefined curve profile
                A_PInv = linalg.pinv(A)                             #Calculate LeftPseudoInverse of predefinedcurve

                LESOutput = dot(A_PInv,CompWindowBuffer)            #Multiply: pinv(A)*Y
                
                #Calculate Phasor Attributes
                Theta_Rad = atan2(LESOutput[1],LESOutput[0])
                Theta_Deg = degrees(Theta_Rad)
                
                Mag = sqrt(((LESOutput[0])**2)+((LESOutput[1])**2)) # ** indicates '^' , to power off.

                

                #Calculating phase angle and frequency change
                #Just for first cycle
                if (x==0): 
                    Theta_Change=0
                    Prev_Theta_Rad = 0
                    if (i==0):
                            csvwrite.writerow(["Time","Magnitude","PhaseAngle","EstimatedFreq"])        #Results headers

                #Every other cycle
                else:
                                                
                        Theta_Change=WrapDifference(Theta_Rad,Prev_Theta_Rad)   #Find Theta difference between +-pi limits
                        Freq=50+((Theta_Change)*Fs/(2*pi*k))                    #Calculate new frequency


        csvwrite.writerow([TimeStamp,Mag,Theta_Rad,Freq])                # write results to csv file
        

        print "Time:{:10.2f}\t  Phasor:{:10.4f}<{:10.4f}\t  Freq:{:10.4f}Hz\t".format(TimeStamp,Mag,Theta_Rad,Freq)
        #print "Time:{:10.2f}\t EstimatedFreq:{:10.5f}Hz\t Theta_Change:{:10.5f}Rad".format(TimeStamp,NewFreq,Theta_Change),"Theta",Theta_Rad,", Prev:",Prev_Theta_Rad,"   ",LESOutput[1],LESOutput[0],Mag
        #print "Time: ",TimeStamp,"  EstFreq:  ",Freq
        Prev_Theta_Rad = Theta_Rad                      #Set 'Prev_Rad' variable

    
def Algorithm_UDP(UDPPort,NumHarmonics):

    UDP_IP = "127.0.0.1"

    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((UDP_IP, UDPPort))
    
    csvwrite = csv.writer(open("Phasor_UDP","wb"))

    Samples_Buffer = 257 #Change this to recieve from UDP?      #Size of computational window
    m = (Samples_Buffer-1)/2                                    #Phasor estimated for middle of computational window
    StartOffset = m                                             #Put centre of computational window at start of 2nd cycle
    Freq = int(50)
    k=(Samples_Buffer-1)



    for x in range(0,3000):     ##change this to loop indefinitely? For payload recieved?
        ##Collect XML Items from UDP
        XMLBuffer, addr = sock.recvfrom(3000) # buffer size is 3000 bytes
        root = ET.fromstring(XMLBuffer)
    
        ##Collect Sampling Frequency
        Fs = float(root[3].text)
        
        ##Collect Timestamp
        TimeStamp=(root[1].text)
        ##Collect Payload and convert - set as "CompWindowBuffer"
        CompWindowBuffer = root[7][4].text
        CompWindowBuffer = PayloadConvert(CompWindowBuffer)
        
        CompWindowBuffer.extend([1]) ##TEMPORARY WORKAROUND! <--
        CompWindowBuffer[(Samples_Buffer-1)]=CompWindowBuffer[(Samples_Buffer-2)]


        
        ##Convert to signed integers:
        i=0
        while i<len(CompWindowBuffer):
                if CompWindowBuffer[i] > 32767:
                        CompWindowBuffer[i] = (65536-CompWindowBuffer[i])*-1
                i +=1

        
        for i in range(0,2):
                A = GenPredefinedCurve(m,NumHarmonics,Fs,Freq)      #Generate predefined curve profile
                A_PInv = linalg.pinv(A)                             #Calculate LeftPseudoInverse of predefinedcurve
                
                LESOutput = dot(A_PInv,CompWindowBuffer)            #Multiply: pinv(A)*Y


                #Calculate Phasor Attributes
                Theta_Rad = atan2(LESOutput[1],LESOutput[0])
                Theta_Deg = degrees(Theta_Rad)

                Mag = sqrt(((LESOutput[0])**2)+((LESOutput[1])**2)) # ** indicates '^' , to power off.

                

                #Calculating phase angle and frequency change
                #Just for first cycle
                if (x==0): 
                    Theta_Change = 0
                    Prev_Theta_Rad=0
                    if(i==0):
                            csvwrite.writerow(["Time","Magnitude","PhaseAngle","EstimatedFreq"])        #Results headers
                            

                #Every other cycle
                else:
                    Theta_Change =WrapDifference(Theta_Rad,Prev_Theta_Rad)                            #phase angle change this window                          
                    Freq=50+((Theta_Change)*Fs/(2*pi*k))                                                #Estimate frequency. = Nominal Freq+Change in Freq


        csvwrite.writerow([TimeStamp,Mag,Theta_Rad,Freq])
        #print "Time: ",TimeStamp,"  EstFreq:  ",Freq
        print "Time:",TimeStamp,"\t  Phasor:{:10.4f}<{:10.4f}\t  Freq:{:10.4f}Hz\t".format(Mag,Theta_Rad,Freq)
        
        Prev_Theta_Rad = Theta_Rad


    





