#/bin/bash

# Delete older versions
rm source/benchmarks/*
rm source/models/*

# Create index file
cat > source/benchmarks/index.rst << EOF
.. include:: ../../../benchmarks/README.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cat > source/models/index.rst << EOF
.. include:: ../../../benchmarks/README.md
   :parser: myst_parser.sphinx_

For more details on the models see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

# Loop over all benchmarks
for f in $(find ../benchmarks/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "benchmarks" ] || continue
    echo "Name of generated file " $NAME.rst
    cat > source/benchmarks/$NAME.rst << EOF
.. include:: ../../$f
   :parser: myst_parser.sphinx_

EOF
    # Append to toctree
    cat >> source/benchmarks/index.rst << EOF
   $NAME
EOF
done

# Loop over all models
for f in $(find ../models/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "models" ] || continue
    echo "Name of generated file " $NAME.rst
    cat > source/models/$NAME.rst << EOF
.. include:: ../../$f
   :parser: myst_parser.sphinx_

EOF
    # Append to toctree
    cat >> source/models/index.rst << EOF
   $NAME
EOF
done
