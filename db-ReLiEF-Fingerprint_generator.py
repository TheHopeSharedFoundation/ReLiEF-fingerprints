################################################
# Copyright 2023 Benjamin M. Samudio
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Ben Samudio, May 2023
# Towards Alleviating Suffering
###############################################
import re
import math
import statistics
import os
import matplotlib

path = "/Users/benjaminsamudio/ReLiEF_Fingerprints_DistributionBased_PDE-10/" #<--------------------------------- This should be set to the path of the working directory which contains the surface dot *.wrl files
os.chdir(path)

rgb_distribution_array = []
rgb_match_count = 0
longest_fingerprint_length = 0

############################################################################################### Section 1 of 3: Create the RGB binner
# The RGB values from the *.wrl file are fractions of one.
rgb_distribution_row_count = 0
rgb_temp_array = []

for rgb_red_value in range(0,52): # <-------------The maximum RGB value is 255.  Here the values are binned in increments of 5-RGB units.
        for rgb_green_value in range(0,52):
                for rgb_blue_value in range(0,52):
                        rgb_distribution_row_count += 1
                        rgb_temp_array = []
                        rgb_bin_string = ""
                        rgb_red_fraction = 5 * rgb_red_value / 255
                        rgb_green_fraction = 5 * rgb_green_value / 255
                        rgb_blue_fraction = 5 * rgb_blue_value / 255
                        rgb_red_top = 5 * (rgb_red_value + 1) / 255
                        rgb_green_top = 5 * (rgb_green_value + 1) / 255
                        rgb_blue_top = 5 * (rgb_blue_value + 1) / 255
                        rgb_middle_hex_code = matplotlib.colors.to_hex([rgb_red_fraction, rgb_green_fraction, rgb_blue_fraction])
                        rgb_bin_string = "PC&" + rgb_middle_hex_code
                        rgb_temp_array = [rgb_distribution_row_count,rgb_red_fraction,rgb_green_fraction,rgb_blue_fraction,rgb_bin_string,0,rgb_red_top,rgb_green_top,rgb_blue_top]
 #                      print(rgb_temp_array)
                        rgb_distribution_array.append(rgb_temp_array)

rgb_distribution_array_size = len(rgb_distribution_array)

fingerprints_filename = path + "ReLieF_Fingerprints_DistributionBased.csv"
header_string = "Name" + "," + "Color" + "," + "Fingerprints" + ","
with open(fingerprints_filename,'w') as file_object:
    header_string = header_string + "\n"
    file_object.write(header_string)




for file in os.listdir(): # START OF ITERATIONS
        if file.endswith(".wrl"):
                print(f"Currently processing this file: {file}")
                temp_filename_tuple = ()
                temp_filename_tuple = os.path.splitext(file)
                output_filename_base = temp_filename_tuple[0]
                dot_surface_file = ""
                dot_surface_file = file
                fingerprint_length = 0
                fingerprint_length_with_delimiter = 0
                fingerprint_output_string = ""
                fingerprint_output_string = output_filename_base + "," + "FILL_ME" + ","
                ############################################################################################### Section 2 of 3: Extract the coordinates and colors of each surface dot.  Convert the colors to an alphabetic code.  Collect all values into: surface_all_dot_values.
                # Example string for coordinates from *.wrl file: translation 2.819419 0.448916 -1.255814
                # Example string for colors from *.wrl file: material Material { diffuseColor 0.0000 0.0624 1.0000

                with open(dot_surface_file) as file_object:
                        for wrl_file_row in file_object: #<-------------------------------------------------- Open a *.wrl file.  This file should only cotain spheres (dots), their translations, and their colors
                                surface_dot_coordinates = []
                                surface_dot_colors = []
                                surface_dot_start = []
                                surface_dot_end = []
                                surface_dot_coordinates = re.search(r"translation\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",wrl_file_row) #<---------------------- Extract dot translation coordinates
                                surface_dot_colors = re.search(r"diffuseColor\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",wrl_file_row) #<-------------------------- Extract dot color
                                surface_dot_start = re.search(r"Transform",wrl_file_row) #<--------------------------------------------------------------------------- Extract the indicator of the START of a new dot block in the *.wrl file
                                surface_dot_end = re.search(r"shininess",wrl_file_row) #<----------------------------------------------------------------------------- Extract the indicator of the END of a new dot block in the *.wrl file
                                if surface_dot_start: #<------------------------------------------- Initiate arrays upon finding a new dot block
                                        new_dot_flag = 1
                                        surface_single_dot_values = []
                                if surface_dot_end:
                                        new_dot_flag = 0
                                if surface_dot_coordinates and new_dot_flag == 1: #<-------------------------------------- These are the dot coordinates (translations) from the *.wrl file.  They're only used for calculating distances.
                                        surface_single_dot_values.append(float(surface_dot_coordinates.group(1)))
                                        surface_single_dot_values.append(float(surface_dot_coordinates.group(2)))
                                        surface_single_dot_values.append(float(surface_dot_coordinates.group(3)))
                                if surface_dot_colors and new_dot_flag == 1: #<-------------------------------------------- These are the RGB values (fractions of 255) from the *.wrl file.
                                        first_color_value = 0
                                        second_color_value = 0
                                        third_color_value = 0
                                        first_color_value = float(surface_dot_colors.group(1))
                                        second_color_value = float(surface_dot_colors.group(2))
                                        third_color_value = float(surface_dot_colors.group(3))
                #                       print(first_color_value,second_color_value,third_color_value)
                                        for rgb_bin_count in range(0,rgb_distribution_array_size):
                                                if first_color_value == 1 and second_color_value == 1 and third_color_value == 1:
                                                        rgb_distribution_array[rgb_bin_count][5] += 1
                                                        rgb_match_count += 1
                #                                       print(rgb_match_count, rgb_distribution_array[rgb_bin_count])
                                                        break
                                                elif first_color_value >= rgb_distribution_array[rgb_bin_count][1] and second_color_value >= rgb_distribution_array[rgb_bin_count][2] and third_color_value >= rgb_distribution_array[rgb_bin_count][3] and first_color_value < rgb_distribution_array[rgb_bin_count][6] and second_color_value < rgb_distribution_array[rgb_bin_count][7] and third_color_value < rgb_distribution_array[rgb_bin_count][8]:
                                                        rgb_match_count += 1
                 #                                      print(rgb_match_count,rgb_distribution_array[rgb_bin_count])
                                                        rgb_distribution_array[rgb_bin_count][5] += 1
                                                        break

                ############################################################################################### Section 3 of 3: Create fingerprints
                #                                                                   bin string
                #                                  V----------------------------------------------------------------------------V
                #                                 property         color bin       match instance           bit            space delimiter
                #                                    V                 V                 V              V----------V             V
                # Example output from this section: PC&#ffa000&1 PC&#ffa000&2 PC&#ffa000&3 PC&#ffa000&4 PC&#ffa000&5 PC&#ffa000&6
                # For each color bin in rgb_distribution_array, the number of matches with dot colors (input RGB values) is tallied.  This tally is the "match/bin instance".
                # The letters "PC" in the example output stand for "partial charge".  The text between the "&" symbols is the color hex code representing the color bin.
                # The number after the last "&" symbol is the match instance.  For each color bin and each instance, a fingerprint "bit" is printed.
                # The example above represents a "bit string" of all of the instances (six in total) of matches that occurred in the partial charge property distribition and color bin "#ffa000".
                # All bit strings are concatenated together into a fingerprint representing a particular surface.

                output_bin_string = ""
                output_fingerprint_string = ""
                output_fingerprint_bits = []
                number_bins_matched = 0
                fingerprint_length = 0
                fingerprint_length_with_delimiter = 0

                for rgb_bin_index in range(0,rgb_distribution_array_size):
                #       output_bin_string = "%" + " " #<------------------------------ Uncomment this only if delimiters are desired between bin strings
                        output_bin_string = ""
                        if rgb_distribution_array[rgb_bin_index][5] != 0: #<---------- Ensure that at least one match was made between a dot color (input RGB values) and a color bin (rgb_distribution_array)
                                rgb_bin_instance = 0 #<---------- This is the match instance count initiation
                                for rgb_bit_index in range(0,rgb_distribution_array[rgb_bin_index][5]):
                                        fingerprint_length += 1
                                        rgb_bin_instance += 1 #<---------- This is the match instance count
                                        output_bin_string = output_bin_string + rgb_distribution_array[rgb_bin_index][4] + "&" + str(rgb_bin_instance) + " "
                                output_bin_string = output_bin_string + "% % % % " #<------------------------------ This adds a buffer to the end of each bin string.
                #               print(output_bin_string) # <---------------------------- Prints bin strings, one per row.
                                output_fingerprint_bits.append(output_bin_string)
                        rgb_distribution_array[rgb_bin_index][5] = 0 #<------------------ Reset the match count in rgb_distribution_array
                number_bins_matched = len(output_fingerprint_bits)
                output_fingerprint_string = "".join(output_fingerprint_bits)
                fingerprint_length_with_delimiter = fingerprint_length + number_bins_matched * 4#<-------------------------- The number_bins_matched is multiplied by 4 to account for buffers
                if fingerprint_length_with_delimiter > longest_fingerprint_length:
                        longest_fingerprint_length = fingerprint_length_with_delimiter
                print(output_fingerprint_string)
                fingerprint_output_string = fingerprint_output_string + output_fingerprint_string
                with open(fingerprints_filename, 'a') as file_object:
                        file_object.write(fingerprint_output_string)
                        file_object.write("\n")
                print(f"The number of color bins matched: {number_bins_matched}")
                print(f"The fingerprint length is: {fingerprint_length_with_delimiter}")
                print(f"The longest fingerprint length is: {longest_fingerprint_length}")
