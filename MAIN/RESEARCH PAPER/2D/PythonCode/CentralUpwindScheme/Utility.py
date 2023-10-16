"""
    Utility module containing useful functions 
"""

def minmod(*nums: float) -> float:
    """
        function implements the minmod limiter
    """
    if min(nums) > 0:
        return min(nums)
    elif max(nums) < 0:
        return max(nums)
    else:
        return 0;

def merge_dicts(dict1: dict, dict2: dict) -> dict:
    """
        Function to merge two dictionaries into one
    """

    return {**dict1, **dict2}


