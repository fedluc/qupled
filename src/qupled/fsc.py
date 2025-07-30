import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

def coth(x):
    return np.cosh(x) / np.sinh(x)

def get_user_input():
    """Get all required input parameters interactively from the user."""
    print("\n" + "="*50)
    print("Finite Size Correction Calculator")
    print("="*50 + "\n")
    
    # Get particle numbers
    n_input = input("Enter particle numbers (comma separated, e.g. 10,100,1000): ")
    N_values = sorted([int(n.strip()) for n in n_input.split(',')])
    
    # Get rs values
    rs_input = input("Enter rs values (comma separated, e.g. 1,5,10): ")
    rs_values = sorted([float(rs.strip()) for rs in rs_input.split(',')])
    
    # Get theta values
    theta_input = input("Enter theta values (comma separated, e.g. 0.5,1,2): ")
    theta_values = sorted([float(theta.strip()) for theta in theta_input.split(',')])
    
    # Fixed G value (maximum value for discrete sum argument and continuous integration)
    G = 50
    
    # Get output directory
    output_dir = input("Enter output directory (press Enter for current directory): ").strip()
    if not output_dir:
        output_dir = "."
    
    return N_values, rs_values, theta_values, G, output_dir

def find_closest_parameters(rs, theta):
    """Find closest existing rs/theta combination in database."""
    all_runs = DataBase.inspect_runs()
    closest_rs = None
    closest_theta = None
    min_distance = float('inf')
    
    for run_info in all_runs:
        db_rs = run_info.get('coupling', None)
        db_theta = run_info.get('degeneracy', None)
        if db_rs is not None and db_theta is not None:
            distance = np.sqrt((db_rs - rs)**2 + (db_theta - theta)**2)
            if distance < min_distance:
                min_distance = distance
                closest_rs = db_rs
                closest_theta = db_theta
    
    return closest_rs, closest_theta

def save_results(all_results, rs, theta, output_dir):
    """Save all results for one rs,theta combination to a single file."""
    filename = f"FSC_results_rs{rs}_theta{theta}.dat"
    filepath = os.path.join(output_dir, filename)
    
    with open(filepath, 'w') as f:
        f.write(f"# Finite Size Correction Results (G=50)\n")
        f.write(f"# rs = {rs}, theta = {theta}\n")
        f.write(f"# {'='*50}\n")
        f.write("# N\tuint_correction\tfxc_correction\n")
        
        for N, results in sorted(all_results.items()):
            f.write(f"{N}\t{results['uint_correction']:.8e}\t{results['fxc_correction']:.8e}\n")
    
    print(f"Results saved to {filepath}")

def compute_fxc_correction(N, theta, G, target_rs, rs_grid=None):
    """
    Compute FSC for both the interaction energy and the exchange correlation
    energy. The sum over the discrete reciprocal vectors is computed only for 
    positive values of l1,l2,l3 due to the negative values being symmetric, i.e.
    we care only about the +++ octant. Then checks are made for how many zeros 
    are present in each combination. In the case that only one of 
    the vectors is non-zero (+00) then the resulting value is multiplied by 2, 
    if two are non-zero (++0) then it's multiplied by 4, and if all three are 
    non-zero (+++) then it's multiplied by 8. Additionally, for further speedup
    we take advantage of NumPy vectorized operations by precomputing the norms 
    for each value of l1, i.e. each 2D slice of the l2,l3 values. 
    """

    # Create rs grid including 0 for f_xc integration, ensuring target_rs is included
    if rs_grid is None:
        # Generate base grid with step 0.1 up to the next multiple of 0.1 above target_rs
        base = np.arange(0.0, target_rs + 0.1, 0.1)
        # Filter points <= target_rs
        base = base[base <= target_rs]
        # Ensure target_rs is included (handles non-multiples of 0.1)
        if base.size == 0 or base[-1] != target_rs:
            base = np.append(base, target_rs)
        rs_grid = base
    
    # Get STLS data for all rs values
    stls_data_dict = {}
    for rs in rs_grid:
        stls_data = find_existing_result(rs, theta)
        stls_data_dict[rs] = stls_data
    
    # Compute rs*uint for all rs values
    rs_uint = []
    
    for rs in rs_grid:
        lambda_val = (4 / (9 * np.pi)) ** (1/3)
        
        # Get S(q) interpolation
        S_interp = interp1d(stls_data_dict[rs]['wvg'],
                           stls_data_dict[rs]['ssf'],
                           kind='cubic',
                           bounds_error=False,
                           fill_value=(1.0, stls_data_dict[rs]['ssf'][-1]))
        
        # =========================================
        # Compute continuous term 
        # =========================================
        def integrand(q):
            return (S_interp(q) - 1)
        
        integration_cutoff = min(G, stls_data_dict[rs]['wvg'][-1])
        continuous = quad(integrand, 0, integration_cutoff)[0]
        continuous *= 1 / (np.pi * lambda_val)
        
        # =========================================
        # Compute discrete term
        # =========================================
        max_l = int(np.ceil((G * (N ** (1/3)) / ((8 * np.pi / 3) ** (1/3))) / np.sqrt(3)))
        prefactor = (2 / (3 * np.pi)) * ((3 / (8 * np.pi)) ** (2/3)) * (1 / (lambda_val)) * (1 / (N ** (1/3)))
        
        # Discrete summation
        total_sum = 0.0
        l2_vals = np.arange(0, max_l + 1)
        l3_vals = np.arange(0, max_l + 1)
        L2, L3 = np.meshgrid(l2_vals, l3_vals, indexing='ij')
        arg_S_const = ((8 * np.pi / 3) ** (1/3)) * (1 / (N ** (1/3)))
        
        # Create 2D vector slices for each l1
        for l1 in range(0, max_l + 1):
            l_norm_grid = np.sqrt(l1**2 + L2**2 + L3**2)
            non_zero_count_grid = (np.full_like(L2, l1 > 0).astype(int) + 
                                 (L2 > 0).astype(int) + 
                                 (L3 > 0).astype(int))
            
            # Find the multiplicity depending on number of zeros 
            multiplicity_grid = np.zeros_like(l_norm_grid, dtype=np.float64)
            multiplicity_grid[non_zero_count_grid == 3] = 8
            multiplicity_grid[non_zero_count_grid == 2] = 4
            multiplicity_grid[non_zero_count_grid == 1] = 2
            
            arg_S_grid = arg_S_const * l_norm_grid
            S_val_grid = S_interp(arg_S_grid)
            
            # Retrieve the full reciprocal space sum value
            with np.errstate(divide='ignore', invalid='ignore'):
                term_grid = (S_val_grid - 1) / (l_norm_grid**2) * multiplicity_grid
                term_grid[np.isinf(term_grid) | np.isnan(term_grid)] = 0
            
            if l1 == 0:
                term_grid[0, 0] = 0.0
            
            total_sum += np.sum(term_grid)
        
        discrete = prefactor * total_sum
        
        # ============================================
        # Compute madelung term
        # ============================================
        madelung = -2.837297 * ((3/(4 * np.pi))**(1/3)) * (N**(-1/3))
        
        # Compute rs*uint
        rs_uint.append(continuous - discrete - madelung/2)

    # Compute the uint FSC
    uint_correction = continuous/target_rs - discrete/target_rs - madelung/2/target_rs
    print(continuous/target_rs)
    print(discrete/target_rs)
    print(madelung/2/target_rs)
    # Integrate rs*uint from 0 to target_rs
    uint_interp = interp1d(rs_grid, rs_uint, kind='cubic', fill_value="extrapolate")
    integral = quad(uint_interp, 0, target_rs)[0]
    # Compute the f_xc FSC
    fxc_correction = integral / (target_rs**2)
    
    return {
        'fxc_correction': fxc_correction,
        'uint_correction': uint_correction,  
        'rs_grid': rs_grid,
        'rs_uint': rs_uint
    }


def run_computations(N_values, rs_values, theta_values, G, output_dir):
    """Run computations and save results for all parameter combinations."""
    lambda_val = (4 / (9 * np.pi)) ** (1/3)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    for rs in rs_values:
        for theta in theta_values:
            print(f"\nProcessing rs = {rs}, theta = {theta}")
            
            # Check if this combination exists
            stls_data = find_existing_result(rs, theta)
            if stls_data is None:
                closest_rs, closest_theta = find_closest_parameters(rs, theta)
                print(f"\nERROR: No STLS data found for rs={rs}, theta={theta}")
                print(f"Closest available combination: rs={closest_rs}, theta={closest_theta}")
                continue
            
            # Store all results for this rs,theta combination
            all_results = {}
            
            # Precompute BCDC values
            BCDC = {N: 1/(4*N)*np.sqrt(3/rs**3)*coth(2/theta*np.sqrt((lambda_val*rs)/(3*np.pi))) 
                   for N in N_values}
            
            # Compute for all N values
            uint_corrections = []
            fxc_corrections = []
            
            for N in N_values:
                print(f"  Computing N = {N}")
                results = compute_fxc_correction(N, theta, G, rs)
                all_results[N] = results
                uint_corrections.append(results['uint_correction'])
                fxc_corrections.append(results['fxc_correction'])
            
            # Save all results for this rs,theta combination
            save_results(all_results, rs, theta, output_dir)
            
            # Generate plot
            plt.figure(figsize=(12, 6))
            
            # Plot uint correction
            plt.subplot(1, 2, 1)
            plt.loglog(1/np.array(N_values), np.abs(uint_corrections), 'o-', 
                     label='FSC', markersize=6)
            plt.loglog(1/np.array(N_values), [abs(BCDC[N]) for N in N_values], 'k--',
                     label='BCDC', linewidth=2)
            plt.title(f'Interaction Energy (rs={rs}, θ={theta})')
            plt.xlabel('1/N')
            plt.ylabel('|FSC|')
            plt.grid(True)
            plt.legend()
            
            # Plot fxc correction
            plt.subplot(1, 2, 2)
            plt.loglog(1/np.array(N_values), np.abs(fxc_corrections), 'o-', 
                     label='FSC', markersize=6)
            plt.title(f'Exchange-Correlation (rs={rs}, θ={theta})')
            plt.xlabel('1/N')
            plt.ylabel('|FSC|')
            plt.grid(True)
            plt.legend()
            
            plt.tight_layout()
            plot_path = os.path.join(output_dir, f'corrections_rs{rs}_theta{theta}.png')
            plt.savefig(plot_path)
            plt.close()
            print(f"  Plot saved to {plot_path}")

def find_existing_result(rs, theta):
    """Check if results for given rs and theta already exist in the database."""
    all_runs = DataBase.inspect_runs()
    
    for run_info in all_runs:
        # Check if parameters match
        if (np.isclose(run_info.get('coupling', None), rs) and 
           np.isclose(run_info.get('degeneracy', None), theta)):
            return DataBase.read_results(run_info['id'])
    
    return None

if __name__ == "__main__":
    # Get all inputs from user
    N_values, rs_values, theta_values, G, output_dir = get_user_input()
    
    # Track if any computations failed
    all_successful = True
    
    # Run computations with user-provided parameters
    for rs in rs_values:
        for theta in theta_values:
            # Check if this combination exists before running computations
            stls_data = find_existing_result(rs, theta)
            if stls_data is None:
                closest_rs, closest_theta = find_closest_parameters(rs, theta)
                print(f"\nERROR: No STLS data found for rs={rs}, theta={theta}")
                print(f"Closest available combination: rs={closest_rs}, theta={closest_theta}")
                all_successful = False
                continue
            
            # Parameters exist run computations
            try:
                run_computations(N_values, [rs], [theta], G, output_dir)
            except Exception as e:
                print(f"\nERROR: Computation failed for rs={rs}, theta={theta}")
                print(f"Error: {str(e)}")
                all_successful = False
    
    print("\n" + "="*50)
    if all_successful:
        print("All computations completed successfully!")
        print(f"Results saved in: {os.path.abspath(output_dir)}")
    else:
        print("Please re-run script.")
    
    print("="*50 + "\n")