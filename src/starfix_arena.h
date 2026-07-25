#ifndef STARFIX_ARENA_H
#define STARFIX_ARENA_H

#ifndef _Alignof
#ifdef __cplusplus
#define _Alignof alignof
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
/* C11 supports _Alignof */
#else
#define _Alignof(type)   \
    offsetof(            \
        struct {         \
            char c;      \
            type member; \
        },               \
        member)
#endif
#endif

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "starfix_status.h"

typedef struct {
    uint8_t* beg;
    uint8_t* end;
} starfix_arena_t;

static inline void* starfix_arena_alloc(starfix_arena_t* a, ptrdiff_t size, ptrdiff_t align,
                                        ptrdiff_t count) {
    if (!a || !a->beg || !a->end || size <= 0 || count < 0) return NULL;
    if (align < 1) align = 1;
    ptrdiff_t padding = (ptrdiff_t)(-(uintptr_t)a->beg & (uintptr_t)(align - 1));
    ptrdiff_t available = (ptrdiff_t)(a->end - a->beg) - padding;
    if (available < 0 || count > available / size) return NULL;
    void* p = a->beg + padding;
    a->beg += padding + count * size;
    return memset(p, 0, (size_t)(count * size));
}

#define starfix_alloc_array(a, type, count) \
    (type*)starfix_arena_alloc((a), sizeof(type), _Alignof(type), (count))

#endif
