
                Case ('KONDO')
                    call Ham_Alloc_Kondo
                Case ('HUBBARD')
                    call Ham_Alloc_Hubbard
                Case ('HUBBARD_PLAIN_VANILLA')
                    call Ham_Alloc_Hubbard_Plain_Vanilla
                Case ('TV')
                    call Ham_Alloc_tV
                Case ('LRC')
                    call Ham_Alloc_LRC
                Case ('Z2_MATTER')
                    call Ham_Alloc_Z2_Matter
                Case ('SPIN_PEIERLS')
                    call Ham_Alloc_Spin_Peierls