import subprocess

def run_orthofinder(input_dir, pixi_env="pixi", additional_args=None, debug=False):
    """
    Run OrthoFinder on a given directory using pixi for environment setup.

    Args:
        input_dir (str): Path to the directory containing input files for OrthoFinder.
        pixi_env (str): The pixi environment command (default is 'pixi').
        additional_args (list): Additional arguments to pass to the OrthoFinder command.
        debug (bool): Whether to print debug information.

    Returns:
        None
    """
    # Ensure additional_args is a list
    if additional_args is None:
        additional_args = []

    # Construct the OrthoFinder command
    command = [pixi_env, "run", "orthofinder", "-f", input_dir] + additional_args

    if debug:
        print(f"Running OrthoFinder with the following command:")
        print(" ".join(command))
        print(f"Input directory: {input_dir}")
        if additional_args:
            print(f"Additional arguments: {additional_args}")

    try:
        # Run the command and capture the output
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Stream the output to the console
        for line in process.stdout:
            print(line, end="")  # Print each line as it is received

        # Wait for the process to complete and capture the return code
        process.wait()

        if process.returncode != 0:
            # If the process failed, print the error output
            error_output = process.stderr.read()
            print(f"\nError: OrthoFinder failed with return code {process.returncode}")
            print(f"Error details:\n{error_output}")
        else:
            if debug:
                print("\nOrthoFinder completed successfully.")

    except FileNotFoundError as e:
        print(f"Error: The command '{command[0]}' was not found. Is pixi installed and in your PATH?")
        if debug:
            print(f"Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred while running OrthoFinder.")
        if debug:
            print(f"Details: {e}")