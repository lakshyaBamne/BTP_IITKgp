# Python file to run the 1D CU code given by Prof. Garg

import numpy as np
import matplotlib.pyplot as plt

def read_output(file_name: str) -> tuple:
    """
        Function to read the output file
    """
    mainlist = []

    with open(file_name, "r") as f:
        data = f.readlines()

    for line in data:
        onelist = line.split(" ")
        temp = []
        for i in onelist:
            if i!='' and i!='\n':
                temp.append(i)
        mainlist.append(temp)
        
    domain = [float(i[0]) for i in mainlist]
    value = [float(i[1]) for i in mainlist]

    return (domain,value)

# def plot(result1: tuple, result2: tuple, result3: tuple) -> None:
#     """
#         Function to plot the results
#     """
#     plt.plot(result1[0], result1[1], color="blue", label="CU")
#     plt.plot(result2[0], result2[1], color="red", linestyle='--',label="CURH") 
#     plt.plot(result3[0], result3[1], color="green", linestyle='-.', linewidth=0.7,label="REFERENCE") 
    
#     plt.legend()
#     plt.savefig("BLW-CU-CURH")
#     plt.show()

def plot(result1: tuple) -> None:
    """
        Function to plot the results
    """
    plt.plot(result1[0], result1[1], color="blue", label="CU")    
    
    plt.legend()
    plt.savefig("BLW-CU-CURH")
    plt.show()

# MCW
# result1 = read_output("out1DMCWrho-CU001")
# result2 = read_output("out1DMCWrho-CURH001")
# result3 = read_output("out1DMCWrho-CU-ref001")

# test
# result1 = read_output("out1Dtest41rho2000")
# result2 = read_output("out1Dtest41rho2001")

# BLW
result1 = read_output("Blast001")
# result2 = read_output("out1DBlastrho-CURH001")
# result3 = read_output("out1DBlastrho-CU-ref001")

# Shu Osher
# result1 = read_output("out1DShuOSrho-CU001")
# result2 = read_output("out1DShuOSrho-CURH001")
# result3 = read_output("out1DShuOsrho-CU-ref001")

# SOD
# result1 = read_output("out1DSODrho-CU001")
# result2 = read_output("out1DSODrho-CURH001")
# result3 = read_output("out1DSODrho-CU-ref001")

# SPP
# result1 = read_output("out1DSPPrho001")
# result2 = read_output("out1DSPPrhoR001")

# Lax
# result1 = read_output("out1DLaxrho-CU001")
# result2 = read_output("out1DLaxrho-CURH001")
# result3 = read_output("out1DLaxrho-CU-ref001")

# plot(result1, result2, result3)

plot(result1)
