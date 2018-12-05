import os

def linear_weight_function(start_time, end_time, current_time):
        """Returns the normalized linear weight for the edge. 
        Keeps decreasing as the current_time moves away from end_time and goes closer to the start_time
        Added 1 in the numerator to avoid weight to go to 0.0
        """
        return (float(current_time - start_time) + 1.0) / float(end_time - start_time)

def get_graph_start_end(in_path):
    """Can be used with the 'create_subgraph_file' function to create a subgraph
    from a specific time slice.
    """
    start_time, end_time = None, None

    with open(read_path, 'r') as edge_list:
        for line in edge_list:
            # Parse the line
            timestamp = int(line.split(' ')[-1])

            # Update start and end time
            if not start_time or start_time > timestamp:
                start_time = timestamp
            if not end_time or end_time < timestamp:
                end_time = timestamp

    return start_time, end_time

def create_subgraph_file(in_path, out_path, start_time, end_time):
    """Example values for creating a 1 year subgraph from the third year:
        in_path = '../data/stackoverflow_full.txt'
        out_path = '../data/stackoverflow_third_year.txt'
        start_time = 1279775877
        end_time = 1310879877
    """
    with open(out_path, 'wb') as output:
        with open(in_path, 'r') as edge_list:
            for line in edge_list:
                # Parse the line
                timestamp = int(line.split(' ')[-1])

                # Copy the line to the output file
                if timestamp >= start_time and timestamp < end_time:
                    output.write(line)
