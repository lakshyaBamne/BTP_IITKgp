"""
    Module to Initialize the variables for the Riemann Problems
"""

def get_initial_conditions(mode: str, problem: str) -> dict:
    """
        Function returns the initial conditions for the Riemann problem
    """

    data = {}

    if problem == "MCW" and mode == "CU":
        data["domainX"] = [-0.2, 0.2]
        data["domainY"] = [0, 0.8]
        data["dx"] = float(1/250)
        data["dy"] = float(1/250)
        data["time"] = [0, 2]
        data["t"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    elif problem == "MCW" and mode == "REF":
        data["domainX"] = [-0.2, 0.2]
        data["domainY"] = [0, 0.8]
        data["dx"] = float(1/1000)
        data["dy"] = float(1/1000)
        data["time"] = [0, 2]
        data["t"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    else:
        print(f"---ERROR--- Please enter correct Problem/Mode ---")


    return data