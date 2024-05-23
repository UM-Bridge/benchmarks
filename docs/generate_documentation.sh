#!/bin/bash
set -eu -o pipefail

mkdir docs_output
cp make.bat Makefile requirements.txt docs_output
mkdir docs_output/markdown
cp -r source docs_output/source
mkdir docs_output/source/forward-benchmarks
mkdir docs_output/source/inverse-benchmarks
mkdir docs_output/source/models
mkdir docs_output/source/umbridge
cp ../README.md docs_output/markdown/main.md
cp ../.readthedocs.yaml docs_output

# Create index files
cp ../benchmarks/README.md docs_output/markdown/forward-benchmarks.md
cat > docs_output/source/forward-benchmarks/index.rst << EOF
==========================
Propagation benchmarks
==========================

.. include:: ../../markdown/forward-benchmarks.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cp ../benchmarks/README.md docs_output/markdown/inverse-benchmarks.md
cat > docs_output/source/inverse-benchmarks/index.rst << EOF
=======================
Inverse benchmarks
=======================

.. include:: ../../markdown/inverse-benchmarks.md
   :parser: myst_parser.sphinx_

For more details on the benchmarks see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cp ../models/README.md docs_output/markdown/models.md
cat > docs_output/source/models/index.rst << EOF
================
Models
================


.. include:: ../../markdown/models.md
   :parser: myst_parser.sphinx_

For more details on the models see their individual documentation:

.. toctree::
   :maxdepth: 3

EOF

cat > docs_output/source/umbridge/index.rst << EOF

.. toctree::
   :maxdepth: 3

   clients
   servers
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

# Loop over UM-Bridge Documentation
cp ../../umbridge/CONTRIBUTING.md docs_output/source/umbridge/contributing.md
for f in $(find ../../umbridge/ -name 'README.md'); do
    NAME=`echo $f | xargs dirname | xargs basename`
    [ "$NAME" != "umbridge"  ] || continue
    cp $f docs_output/source/umbridge/$NAME.md
    [ "$NAME" == "clients" ] || continue
    [ "$NAME" == "servers" ] || continue
    # Append to toctree
    cat >> docs_output/source/umbridge/index.rst << EOF
   $NAME
EOF
done

