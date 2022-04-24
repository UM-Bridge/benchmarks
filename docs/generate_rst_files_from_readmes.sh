#/bin/bash

# Loop over all benchmarks
for f in $(find ../benchmarks/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ ! -z "$NAME" ] || NAME=index
    echo "Name of generated file " $NAME.rst
    cat > source/benchmarks/$NAME.rst << EOF
.. include::../$f
:parser: myst_parser.sphinx_
EOF
done

# Loop over all models
for f in $(find ../models/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ ! -z "$NAME" ] || NAME=index
    echo "Name of generated file " $NAME.rst
    cat > source/models/$NAME.rst << EOF
.. include:: ../$f
    :parser: myst_parser.sphinx_
EOF
done
