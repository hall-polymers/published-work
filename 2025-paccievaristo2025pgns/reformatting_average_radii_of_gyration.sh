#!/opt/homebrew/Cellar/bash/5.2.37/bin/bash

awk '
    BEGIN {
        # Reformatting of comments.
        getline first_line
        print first_line
        getline
        split($0, second_line)
        getline
        split($0, third_line)
        printf "# %s %s\n", second_line[2], third_line[3]

        # Reformatting of data.
        while (getline > 0) {
            timestep = $1
            number_rows = $2
            printf "%s", timestep
            for (index_row = 1; index_row <= number_rows; index_row++) {
                getline
                row_value = $NF
                printf " %s", row_value
            }
            print ""
        }
    }
'
