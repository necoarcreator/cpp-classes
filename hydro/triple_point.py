# -*- coding: utf-8 -*-
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob
import re

# Visualization settings
PLOT_SETTINGS = {
    'figure': {
        'size': (10, 8),
        'dpi': 100
    }
}

def read_data(file_path):
    """Read and parse data file"""
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        raw_data = [line.strip().split(',') for line in f if line.strip()]
    
    data = np.array([[float(x) for x in row] for row in raw_data])
    
    # Find point with maximum p
    max_p_idx = np.argmax(data[2])  # Column 2 is p
    point = {
        'x': data[max_p_idx, 0],
        'y': data[max_p_idx, 1],
        'p': data[max_p_idx, 2],
        'vx': data[max_p_idx, 3],
        'vy': data[max_p_idx, 4]
    }
    
    return point

def extract_time_value(filename):
    """Extract time value from filename like 'out_0.003357_.csv'"""
    match = re.search(r'out_([0-9.]+)_\.csv', filename)
    if match:
        return float(match.group(1))
    return 0.0  # Default value if pattern not found

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # Check arguments
    if len(sys.argv) != 3:
        if rank == 0:
            print("Usage: mpirun -n <procs> python pressure_trajectory.py <input_dir> <output_dir>")
        return
    
    input_dir, output_dir = sys.argv[1], sys.argv[2]
    
    # Create output directory
    if rank == 0:
        os.makedirs(output_dir, exist_ok=True)
    comm.Barrier()
    
    # Get file list
    if rank == 0:
        files = sorted(glob.glob(os.path.join(input_dir, 'out_*.csv')), 
                      key=lambda x: float(re.search(r'out_([0-9.]+)_\.csv', x).group(1)))
        if not files:
            print(f"No CSV files found in {input_dir}")
            return
    else:
        files = None
    
    files = comm.bcast(files, root=0)
    
    # Distribute files
    n_files = len(files)
    chunk = n_files // size
    remainder = n_files % size
    start = rank * chunk + min(rank, remainder)
    end = start + chunk + (1 if rank < remainder else 0)
    local_files = files[start:end]
    
    # Collect trajectory points
    local_trajectory = []
    
    for file_path in local_files:
        try:
            point = read_data(file_path)
            # Extract time value from filename
            time_value = extract_time_value(os.path.basename(file_path))
            point['time'] = time_value
            local_trajectory.append(point)
        except Exception as e:
            print(f"[Rank {rank}] Error processing {file_path}: {str(e)}")
    
    # Gather all trajectory points on rank 0
    all_trajectory = comm.gather(local_trajectory, root=0)
    
    if rank == 0:
        # Combine and sort points by time
        trajectory = [point for sublist in all_trajectory for point in sublist]
        trajectory.sort(key=lambda x: x['time'])
        
        # Create plot
        plt.figure(figsize=PLOT_SETTINGS['figure']['size'], 
                 dpi=PLOT_SETTINGS['figure']['dpi'])
        
        # Extract coordinates
        x_coords = [p['x'] for p in trajectory]
        y_coords = [p['y'] for p in trajectory]
        
        # Plot trajectory
        plt.plot(x_coords, y_coords, 'b-', linewidth=2, label='Trajectory')
        scatter = plt.scatter(x_coords, y_coords, c=[p['p'] for p in trajectory], 
                   cmap='viridis', s=50, label='Pressure points')
        
        # Add direction arrows (only if there's movement)
        if len(trajectory) > 1:
            step = max(1, len(trajectory)//10)
            for i in range(0, len(trajectory)-1, step):
                dx = trajectory[i+1]['x'] - trajectory[i]['x']
                dy = trajectory[i+1]['y'] - trajectory[i]['y']
                # Only draw arrow if there's actual movement
                if dx != 0 or dy != 0:
                    plt.arrow(trajectory[i]['x'], trajectory[i]['y'], 
                             dx*0.9, dy*0.9, 
                             head_width=0.05, head_length=0.1, 
                             fc='r', ec='r')
        
        # Configure plot
        plt.colorbar(scatter, label='Pressure')
        plt.title('Trajectory of Maximum Pressure Point')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.grid(True)
        plt.legend()
        
        # Save result
        output_path = os.path.join(output_dir, 'pressure_trajectory.png')
        plt.savefig(output_path, dpi=PLOT_SETTINGS['figure']['dpi'], 
                   bbox_inches='tight')
        plt.close()
        
        print(f"Trajectory plot saved to {output_path}")

if __name__ == "__main__":
    main()