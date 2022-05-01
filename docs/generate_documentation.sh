#/bin/bash
set -eu -o pipefail

mkdir docs_output
cp make.bat Makefile requirements.txt docs_output
mkdir docs_output/markdown
cp -r source docs_output/source
cp ../README.md docs_output/markdown/main.md

# Create index file
cp ../benchmarks/README.md docs_output/markdown/benchmarks.md
cat > source/benchmarks/index.rst << EOF
.. include:: ../markdown/benchmarks.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cp ../models/README.md docs_output/markdown/models.md
cat > source/models/index.rst << EOF
.. include:: ../markdown/models.md
   :parser: myst_parser.sphinx_

For more details on the models see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

# Loop over all benchmarks
for f in $(find ../benchmarks/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "benchmarks" ] || continue
    cp $f source/benchmarks/$NAME.md
    # Append to toctree
    cat >> source/benchmarks/index.rst << EOF
   $NAME
EOF
done

# Loop over all models
for f in $(find ../models/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "models" ] || continue
    cp $f source/models/$NAME.md
    # Append to toctree
    cat >> source/models/index.rst << EOF
   $NAME
EOF
done
