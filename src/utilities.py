import os

def linear_weight_function(start_time, end_time, current_time):
        """Returns the normalized linear weight for the edge. 
        Keeps decreasing as the current_time moves away from end_time and goes closer to the start_time
        Added 1 in the numerator to avoid weight to go to 0.0
        """
        return (float(current_time - start_time) + 1.0) / float(end_time - start_time)