# 
# This script executes script GF_standalone (compiles and execute GF) and then compare with base results.
#

# MONAN
#
# Author: Denis Eiras
#
# E-mail: denis.eiras@inpe.br  
#
# Date: 2023/04/13
#
# Version: TODO -> USAR VERSÃO DA RELEASE, QUANDO FOR GERADA
#
#
# Full description
# ================
# - executes compare_diff_and_contourn_all_variables.gs script
# - run binay comparison between grads .gra files
# - run grads diff comparison of variables. Generates figures of diff and contour
# - run scripts
#
#
# Instructions
# ~~~~~~~~~~~~
# 
# - First of all:
#
#   - The base results are reference files must be at ../refs.
#
# - Now you can execute this script:
#
#   $ ./compare_local_modifications.sh
#
#
# History
# =======
# 2023/04/13 - Creation - deniseiras
#
#
# Licence
# =======
#
# <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the  Free  Software  Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it  will be useful, but
# ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
# **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
# GNU General Public License for more details.
#
# You should have received a copy  of the GNU General  Public  License
# along with this program.  If not, see [GNU Public    License](https://www.gnu.org/licenses/gpl-3.0.html).
#

echo -e "\n\nCompiling and executing WSM5"
echo -e "==============================\n"

./get_input_data.sh
mv -f WSM05_dataIn-03600-00001.bin ../datain

./WSM5_standalone.sh
echo -e "\n\nComparing binary diferences"
echo -e "===========================\n"
diff -qr ../dataout/check_WSM05_dataOut-03600-00001.gra ../refs/WSM05_dataOut-ref-03600-00001.gra > diff.txt
if [ -s diff.txt ]; then
  # The file is not-empty.
  echo
  echo -e "\n\nComparing using grads"
  echo -e "=====================\n"
  grads -lc "compare_diff_and_contourn_all_variables.gs ../refs/WSM05_dataOut-ref-03600-00001.ctl ../dataout/check_WSM05_dataOut-03600-00001.ctl"

else
  echo -e "\nCongratulations, no diferences found!!!"
fi
