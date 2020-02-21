#!/usr/bin/env bash

[ -d "$HOME/.ccache" ] && sudo chown -R $USER: $HOME/.ccache
export parent_dir=$(pwd)
mkdir -p $parent_dir/Build/linux/${build_type:-Release}
cd $parent_dir/Build/linux/${build_type:-Release}
cmake $parent_dir -G"${generator:-Unix Makefiles}" -DCMAKE_BUILD_TYPE=${build_type:-Release} -DBUILD_SHARED_LIBS=${shared_libs:-ON} -DBUILD_TESTING=${testing:-OFF} ${CMAKE_EFLAGS}
cmake -j$(nproc) --build .
sudo cmake --build . --target install
cd $parent_dir
$parent_dir/Bin/${build_type:-Release}/SvtAv1EncApp -help
