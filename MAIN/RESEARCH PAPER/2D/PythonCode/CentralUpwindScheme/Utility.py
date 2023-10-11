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




