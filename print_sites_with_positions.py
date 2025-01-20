def print_sites_with_positions(sites_list):
    """
    Print a list of sites with their positions horizontally and vertically aligned.
    Input: a list of tuples, each containing an amino acid/nucleotide and a position.
    Output: None, but prints the amino acids horizontally and the positions vertically aligned.
    """
    amino_acids = [site[0] for site in sites_list]
    positions = [str(site[1]) for site in sites_list]

    # Print amino acids horizontally
    print("".join(amino_acids))

    # Find the maximum number of digits in positions
    max_digits = max(len(pos) for pos in positions)

    # Print positions vertically aligned
    for i in range(max_digits):
        for pos in positions:
            if len(pos) > i:
                print(pos[i], end="")
            else:
                print(" ", end="")
        print()


if __name__ == "__main__":
    # Example usage:
    import random

    x = [(random.choice("ATGC"), i + 1) for i in range(60)]
    print_sites_with_positions(x)
