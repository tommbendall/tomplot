"""assorted routines for post plotting processing"""
import subprocess
import os
import re


def create_animation(results_dir, plot_dir, file_name='animation.gif',
                     delay=20, loop=0):
    u"""
    A routine which makes use of the imageMagick function 'convert' to create a
    gif from a directory of image files.
    Requirements: Requires an installation of ImageMagick on your bash path
    Args:
        results_dir (class: string): file path to the file containg the
                                     images to make the gif from.
        plot_dir (class: string): file path to where the gif is outputted.
        file_name (class: string): name of animation, typicaly end in git
        delay (class: int): delay for animation, measured in milliseconds
        loop (class: int): how many times the animation will loop. 0 = infinite
    u"""
    # check if the input directory exists
    if not os.path.exists(results_dir):
        raise ValueError(f'Error: Input directory, {results_dir} does not exist.')

    if not file_name.endswith('.gif') == True:
        raise ValueError(f'Error: file name, {file_name} lacks.gif.')

    image_files = [f for f in os.listdir(plot_dir) if
                   f.lower().endswith(('.png', '.jpg', '.jpeg'))]

    # function to ensure files get sorted in ascending order
    def extract_number(filename):
        return int(re.search(r'\d+', filename).group())

    sorted_files = sorted(image_files, key=extract_number)
    # configures the convert bash command
    convert_cmd = ['convert',
                   '-delay', str(delay),
                   '-loop', str(loop)]
    # add the input files to the command
    convert_cmd.extend([os.path.join(plot_dir, f) for f in sorted_files])

    # add output file
    convert_cmd.append(f'{results_dir}/{file_name}')

    # use subprocesses to call the command
    subprocess.run(convert_cmd)
