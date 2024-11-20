# List all files in the directory except the ones in keep_list.txt
find . -maxdepth 1 -type f | grep -v -f keep_list.txt > to_delete.txt

# Review the files to be deleted
cat to_delete.txt

# If satisfied, delete the listed files
xargs rm < to_delete.txt

# Optional: Clean up the temporary file
#rm to_delete.txt
