mkdir -p tmp

cd ../data

mkdir -p matrix loom

cp -r old1/outs/filtered_feature_bc_matrix/ matrix/old1
cp -r old2/outs/filtered_feature_bc_matrix/ matrix/old2
cp -r old3/outs/filtered_feature_bc_matrix/ matrix/old3
cp -r young1/outs/filtered_feature_bc_matrix/ matrix/young1
cp -r young2/outs/filtered_feature_bc_matrix/ matrix/young2
cp -r young3/outs/filtered_feature_bc_matrix/ matrix/young3

cp -r old1/velocyto/old1.loom loom/
cp -r old2/velocyto/old2.loom loom/
cp -r old3/velocyto/old3.loom loom/
cp -r young1/velocyto/young1.loom loom/
cp -r young2/velocyto/young2.loom loom/
cp -r young3/velocyto/young3.loom loom/
