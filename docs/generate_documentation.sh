#/bin/bash
set -eu -o pipefail

mkdir docs_output
cp make.bat Makefile requirements.txt docs_output
mkdir docs_output/markdown
cp -r source docs_output/source
mkdir docs_output/source/forward-benchmarks
mkdir docs_output/source/inverse-benchmarks
mkdir docs_output/source/models
cp ../README.md docs_output/markdown/main.md

# Create index file
cp ../benchmarks/README.md docs_output/markdown/forward-benchmarks.md
cat > docs_output/source/forward-benchmarks/index.rst << EOF
.. include:: ../markdown/forward-benchmarks.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cp ../benchmarks/README.md docs_output/markdown/inverse-benchmarks.md
cat > docs_output/source/inverse-benchmarks/index.rst << EOF
.. include:: ../markdown/inverse-benchmarks.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cp ../models/README.md docs_output/markdown/models.md
cat > docs_output/source/models/index.rst << EOF
.. include:: ../markdown/models.md
   :parser: myst_parser.sphinx_

For more details on the models see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

# Loop over all forward benchmarks
for f in $(find ../benchmarks/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "benchmarks" ] || continue
    [[ "$NAME" == *"propagation"* ]] || continue
    cp $f docs_output/source/forward-benchmarks/$NAME.md
    # Append to toctree
    cat >> docs_output/source/forward-benchmarks/index.rst << EOF
   $NAME
EOF
done

# Loop over all inverse benchmarks
for f in $(find ../benchmarks/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "benchmarks" ] || continue
    [[ "$NAME" != *"propagation"* ]] || continue
    cp $f docs_output/source/inverse-benchmarks/$NAME.md
    # Append to toctree
    cat >> docs_output/source/inverse-benchmarks/index.rst << EOF
   $NAME
EOF
done

# Loop over all models
for f in $(find ../models/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "models" ] || continue
    cp $f docs_output/source/models/$NAME.md
    # Append to toctree
    cat >> docs_output/source/models/index.rst << EOF
   $NAME
EOF
done
