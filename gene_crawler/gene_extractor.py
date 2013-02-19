#!/usr/bin/env python

def extract_gene(f, output, start, length):
    """
    Takes an input file, goes to a starting location and extracts the specified number of bases.
    """
    
    line_width = 70
    write_counter = 0
    base_counter = 0
    curr_seg = ""
    for line in f:
        if length == 0:
            break
        line = line.strip()
        if not line:
            continue
        if line[0] == ">":
            break
        
        # Check to see if we have a portion of the segment that is after or 
        # contains the starting base
        segment_start = base_counter
        segment_end = base_counter + len(line)
        base_counter += len(line)
        if segment_end <= start:
            continue

        # Extract all or part of the current segment that should be used.
        if segment_start < start:
            curr_seg = line[start-segment_start:]
        else:
            curr_seg += line
            
        # If we have at least a line's worth of data, write it.
        while length and len(curr_seg) >= line_width:
            amount_to_write = line_width
            if amount_to_write > length:
                amount_to_write = length
            output.write("%s\n" % curr_seg[0:amount_to_write])
            length -= amount_to_write
            curr_seg = curr_seg[amount_to_write:]

    # If after the for loop there is any data left over and we have not
    # hit the write limit, write the data.
    if curr_seg:
        if length:
            amount_to_write = line_width
            if amount_to_write > length:
                amount_to_write = length
            output.write("%s\n" % curr_seg[0:amount_to_write])

