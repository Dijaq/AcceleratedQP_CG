#ifdef IGL_INLINE
#undef IGL_INLINE
#endif

#ifndef IGL_STATIC_LIBRARY
#  define IGL_INLINE inline
#else
#  define IGL_INLINE
#endif
