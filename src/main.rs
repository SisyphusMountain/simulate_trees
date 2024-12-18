use newick_parser::node::{FlatTree, flat_to_node};
use newick_parser::newick::node_to_newick;
use rand::thread_rng;
use rand_distr::{Exp, Distribution};
use rand::Rng;
use rand::seq::index::sample;
use std::env;
use std::process;
use std::fs;
use std::time::Instant;

fn conditional_bd(birth_rate: f64, death_rate: f64, n_extant: usize) -> FlatTree {
    let mut tree = FlatTree{nodes: Vec::new(),
                            root: 100000};
    let mut node_name_counter = 0;
    tree.nodes.reserve(2 * n_extant);
    // Initialize the tree with n_extant nodes
    for _ in 0..n_extant {
        tree.add_node(
            node_name_counter.to_string(),
            None,
            None,
            None,
            Some(0.0),
            0.0,
        );
        node_name_counter += 1;
    }

    let mut currently_alive_indices: Vec<usize> = (0..n_extant).collect();
    let mut n = n_extant;
    let mut current_time = 0.0;
    let mut trng = thread_rng();
    let p_birth = birth_rate / (birth_rate + death_rate);
    let total_rate = birth_rate + death_rate;
    let exp_dist = Exp::new(total_rate).unwrap();
    while n > 0 {
        let waiting_time = exp_dist.sample(&mut trng) / (n as f64);
        current_time += waiting_time;

        let is_birth = trng.gen_bool(p_birth);
        if is_birth && n > 1 {
            // Choose two nodes to fuse
            let i1 = trng.gen_range(0..n);
            let mut i2 = trng.gen_range(0..(n - 1));
            if i2 >= i1 {
                i2 += 1;
            }
            
            let left_child_index = currently_alive_indices[i1];
            let right_child_index = currently_alive_indices[i2];
            // Make sure to remove the higher index first to avoid invalidating the lower index
            // Remove them:
            let (max_pos, min_pos) = if i1 > i2 { (i1, i2) } else { (i2, i1) };
            currently_alive_indices.swap_remove(max_pos);
            currently_alive_indices.swap_remove(min_pos);

            let new_node_name = node_name_counter.to_string();
            tree.add_node(
                new_node_name,
                Some(left_child_index),
                Some(right_child_index),
                None,
                Some(current_time),
                0.0,
            );

            tree[left_child_index].parent = Some(node_name_counter);
            tree[left_child_index].length = current_time - tree[left_child_index].depth.expect("Depth not set");
            tree[right_child_index].parent = Some(node_name_counter);
            tree[right_child_index].length = current_time - tree[right_child_index].depth.expect("Depth not set");

            currently_alive_indices.push(node_name_counter);
            node_name_counter += 1;

            

            n -= 1;
        } 
        else if is_birth && n == 1 {
            // Add the length of the last node
            let last_node_index = currently_alive_indices[0];
            tree[last_node_index].length = current_time - tree[last_node_index].depth.expect("Depth not set");
            n -= 1;
        }
        
        else {
            // Add a new node disconnected from the others
            let new_node_name = node_name_counter.to_string();
            tree.add_node(
                new_node_name,
                None,
                None,
                None,
                Some(current_time),
                0.0,
            );
            currently_alive_indices.push(node_name_counter);
            node_name_counter += 1;
            n += 1;
        }
    }
    tree.root = node_name_counter - 1;
    // Now correct all node depths and lengths
    let max_height: f64 = tree[tree.root].depth.expect("Depth not set") + tree[tree.root].length;
    for node in tree.nodes.iter_mut() {
        node.depth = Some(max_height - node.depth.expect("Depth not set"));
    }
    tree
}

fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();

    // Ensure there are exactly three arguments
    if args.len() != 4 {
        eprintln!("Usage: {} <birth_rate> <death_rate> <n_extant>", args[0]);
        process::exit(1);
    }

    // Parse the arguments
    let birth_rate: f64 = args[1].parse().unwrap_or_else(|_| {
        eprintln!("Error: birth_rate must be a floating-point number");
        process::exit(1);
    });
    let death_rate: f64 = args[2].parse().unwrap_or_else(|_| {
        eprintln!("Error: death_rate must be a floating-point number");
        process::exit(1);
    });
    let n_extant: usize = args[3].parse().unwrap_or_else(|_| {
        eprintln!("Error: n_extant must be a positive integer");
        process::exit(1);
    });

    // Start timing
    let start_time = Instant::now();

    // Call the main tree simulation function
    let tree = conditional_bd(birth_rate, death_rate, n_extant);
    let node = flat_to_node(&tree.nodes, tree.root).expect("Node generation from flat tree failed");
    let mut newick: String = node_to_newick(&node);

    // Add a semicolon at the end of the Newick string
    newick.push(';');
    let elapsed_time = start_time.elapsed();
    // Save the Newick string to a file
    let output_file = "tree.nwk";
    fs::write(output_file, &newick).unwrap_or_else(|_| {
        eprintln!("Error: Failed to write to file {}", output_file);
        process::exit(1);
    });



    // Print results
    println!("Newick string saved to '{}'.", output_file);
    println!("Execution time before writing to disk: {:.2?}", elapsed_time);
}