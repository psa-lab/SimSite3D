aclocal
libtoolize --copy --force --automake
automake -a -c
autoconf
./configure --prefix=`pwd`
