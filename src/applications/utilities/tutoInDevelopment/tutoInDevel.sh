echo "set -e; if [ -z \$1 ]; then debug=no; else debug=\$1; fi; if [ -z \$2 ]; then tutoInDevelopment -debug \$debug; else tutoInDevelopment -debug \$debug -errorMsg \"\${2}\"; fi"
