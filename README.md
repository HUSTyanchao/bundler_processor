# bundler_processor

This a c++ project that process the 3D point in the file bundle.out from
the software bundler(https://www.cs.cornell.edu/~snavely/bundler/)

You can convert the bundle.out into ply file so that it can be visualized
by MeshLab(http://meshlabstuff.blogspot.com/)

Also you can do some statistics about the 3D point cloud of your scene.

CMakeLists.txt is provided so that it can be use both on Windows and Linux.


usage:

./bundler_processor + data_path + bundle_out_file + picture_list  + what_kind_of_process

for example, bundler_processor E:/Dubrovnik6K /bundle/bundle.orig.out list.orig.txt 1

in the data_path E:/Dubrovnik6K, find the bundle folder and bundle.orig.out file in it,

with list.orig.txt, and do FIND_QUERY_TRUE_MATCH process.
