__precompile__()
"""
SU_mod

Implements seismic un*x (SU)  format 
and I/O of SU files

"""
module SU_mod


export SU_hdr, SU_trc, SUPutTrace, SUGetTrace

mutable struct SU_hdr
        tracl::Int32                                       #  1-4    
        tracr::Int32                                       #  5-8  
        fldr::Int32                                        #  9-12 
        tracf::Int32                                       # 13-16 
        ep::Int32                                          # 17-20 
        cdp::Int32                                         # 21-24 
        cdpt::Int32                                        # 25-28 
        trid::Int16                                        # 29-30 
        nvs::Int16                                         # 31-32 
        nhs::Int16                                         # 33-34 
        duse::Int16                                        # 35-36 
        offset::Int32                                      # 37-40 
        gelev::Int32                                       # 41-44 
        selev::Int32                                       # 45-48 
        sdepth::Int32                                      # 49-52 
        gdel::Int32                                        # 53-56 
        sdel::Int32                                        # 57-60 
        swdep::Int32                                       # 61-64 
        gwdep::Int32                                       # 65-68 
        scalel::Int16                                      # 69-70 
        scalco::Int16                                      # 71-72 
        sx::Int32                                          # 73-76 
        sy::Int32                                          # 77-80 
        gx::Int32                                          # 81-84 
        gy::Int32                                          # 85-88 
        counit::Int16                                      # 89-90 
        wevel::Int16                                       # 91-92 
        swevel::Int16                                      # 93-94 
        sut::Int16                                         # 95-96 
        gut::Int16                                         # 97-98 
        sstat::Int16                                       # 99-100 
        gstat::Int16                                       #101-102 
        tstat::Int16                                       #103-104 
        laga::Int16                                        #105-106 
        lagb::Int16                                        #107-108 
        delrt::Int16                                       #109-110 
        muts::Int16                                        #111-112 
        mute::Int16                                        #113-114 
        ns::Int16                                          #115-116 
        dt::Int16                                          #117-118 
        gain::Int16                                        #119-120 
        igc::Int16                                         #121-122 
        igi::Int16                                         #123-124 
        corr::Int16                                        #125-126 
        sfs::Int16                                         #127-128 
        sfe::Int16                                         #129-130 
        slen::Int16                                        #131-132 
        styp::Int16                                        #133-134 
        stas::Int16                                        #135-136 
        stae::Int16                                        #137-138 
        tatyp::Int16                                       #139-140 
        afilf::Int16                                       #141-142 
        afils::Int16                                       #143-144 
        nofilf::Int16                                      #145-146 
        nofils::Int16                                      #147-148 
        lcf::Int16                                         #149-150 
        hcf::Int16                                         #151-152 
        lcs::Int16                                         #153-154 
        hcs::Int16                                         #155-156 
        year::Int16                                        #157-158 
        day::Int16                                         #159-160 
        hour::Int16                                        #161-162 
        minute::Int16                                      #163-164 
        sec::Int16                                         #165-166 
        timbas::Int16                                      #167-168 
        trwf::Int16                                        #169-170 
        grnors::Int16                                      #171-172 
        grnofr::Int16                                      #173-174 
        grnlof::Int16                                      #175-176 
        gaps::Int16                                        #177-178 
        otrav::Int16                                       #179-180 
        d1::Float32                                        #181-184 
        f1::Float32                                        #185-188 
        d2::Float32                                        #189-192 
        f2::Float32                                        #193-196 
        ungpow::Float32                                    #197-200
        unscale::Float32                                   #201-204 
        ntr::Int32                                         #205-208 
        mark::Int16                                        #209-210 
        shortpad::Int16                                    #211-212 
        unassigned::Array{Int16,1}
        #     
        function SU_hdr(
                tracl=Int32(0),                                       #  1-4    
                tracr=Int32(0),                                       #  5-8  
                fldr=Int32(0),                                        #  9-12 
                tracf=Int32(0),                                       # 13-16 
                ep=Int32(0),                                          # 17-20 
                cdp=Int32(0),                                         # 21-24 
                cdpt=Int32(0),                                        # 25-28 
                trid=Int16(0),                                        # 29-30 
                nvs=Int16(0),                                         # 31-32 
                nhs=Int16(0),                                         # 33-34 
                duse=Int16(0),                                        # 35-36 
                offset=Int32(0),                                      # 37-40 
                gelev=Int32(0),                                       # 41-44 
                selev=Int32(0),                                       # 45-48 
                sdepth=Int32(0),                                      # 49-52 
                gdel=Int32(0),                                        # 53-56 
                sdel=Int32(0),                                        # 57-60 
                swdep=Int32(0),                                       # 61-64 
                gwdep=Int32(0),                                       # 65-68 
                scalel=Int16(0),                                      # 69-70 
                scalco=Int16(0),                                      # 71-72 
                sx=Int32(0),                                          # 73-76 
                sy=Int32(0),                                          # 77-80 
                gx=Int32(0),                                          # 81-84 
                gy=Int32(0),                                          # 85-88 
                counit=Int16(0),                                      # 89-90 
                wevel=Int16(0),                                       # 91-92 
                swevel=Int16(0),                                      # 93-94 
                sut=Int16(0),                                         # 95-96 
                gut=Int16(0),                                         # 97-98 
                sstat=Int16(0),                                       # 99-100 
                gstat=Int16(0),                                       #101-102 
                tstat=Int16(0),                                       #103-104 
                laga=Int16(0),                                        #105-106 
                lagb=Int16(0),                                        #107-108 
                delrt=Int16(0),                                       #109-110 
                muts=Int16(0),                                        #111-112 
                mute=Int16(0),                                        #113-114 
                ns=Int16(0),                                          #115-116 
                dt=Int16(0),                                          #117-118 
                gain=Int16(0),                                        #119-120 
                igc=Int16(0),                                         #121-122 
                igi=Int16(0),                                         #123-124 
                corr=Int16(0),                                        #125-126 
                sfs=Int16(0),                                         #127-128 
                sfe=Int16(0),                                         #129-130 
                slen=Int16(0),                                        #131-132 
                styp=Int16(0),                                        #133-134 
                stas=Int16(0),                                        #135-136 
                stae=Int16(0),                                        #137-138 
                tatyp=Int16(0),                                       #139-140 
                afilf=Int16(0),                                       #141-142 
                afils=Int16(0),                                       #143-144 
                nofilf=Int16(0),                                      #145-146 
                nofils=Int16(0),                                      #147-148 
                lcf=Int16(0),                                         #149-150 
                hcf=Int16(0),                                         #151-152 
                lcs=Int16(0),                                         #153-154 
                hcs=Int16(0),                                         #155-156 
                year=Int16(0),                                        #157-158 
                day=Int16(0),                                         #159-160 
                hour=Int16(0),                                        #161-162 
                minute=Int16(0),                                      #163-164 
                sec=Int16(0),                                         #165-166 
                timbas=Int16(0),                                      #167-168 
                trwf=Int16(0),                                        #169-170 
                grnors=Int16(0),                                      #171-172 
                grnofr=Int16(0),                                      #173-174 
                grnlof=Int16(0),                                      #175-176 
                gaps=Int16(0),                                        #177-178 
                otrav=Int16(0),                                       #179-180 
                d1=Float32(0.0),                                      #181-184 
                f1=Float32(0.0),                                      #185-188 
                d2=Float32(0.0),                                      #189-192 
                f2=Float32(0.0),                                      #193-196 
                ungpow=Float32(0.0),                                  #197-200
                unscale=Float32(0.0),                                 #201-204 
                ntr=Int32(0),                                         #205-208 
                mark=Int16(0),                                        #209-210 
                shortpad=Int16(0),                                    #211-212 
                unassigned=zeros(Int16,14))                           #213-240
                new(
                    tracl,                                       #  1-4    
                    tracr,                                       #  5-8  
                    fldr,                                        #  9-12 
                    tracf,                                       # 13-16 
                    ep,                                          # 17-20 
                    cdp,                                         # 21-24 
                    cdpt,                                        # 25-28 
                    trid,                                        # 29-30 
                    nvs,                                         # 31-32 
                    nhs,                                         # 33-34 
                    duse,                                        # 35-36 
                    offset,                                      # 37-40 
                    gelev,                                       # 41-44 
                    selev,                                       # 45-48 
                    sdepth,                                      # 49-52 
                    gdel,                                        # 53-56 
                    sdel,                                        # 57-60 
                    swdep,                                       # 61-64 
                    gwdep,                                       # 65-68 
                    scalel,                                      # 69-70 
                    scalco,                                      # 71-72 
                    sx,                                          # 73-76 
                    sy,                                          # 77-80 
                    gx,                                          # 81-84 
                    gy,                                          # 85-88 
                    counit,                                      # 89-90 
                    wevel,                                       # 91-92 
                    swevel,                                      # 93-94 
                    sut,                                         # 95-96 
                    gut,                                         # 97-98 
                    sstat,                                       # 99-100 
                    gstat,                                       #101-102 
                    tstat,                                       #103-104 
                    laga,                                        #105-106 
                    lagb,                                        #107-108 
                    delrt,                                       #109-110 
                    muts,                                        #111-112 
                    mute,                                        #113-114 
                    ns,                                          #115-116 
                    dt,                                          #117-118 
                    gain,                                        #119-120 
                    igc,                                         #121-122 
                    igi,                                         #123-124 
                    corr,                                        #125-126 
                    sfs,                                         #127-128 
                    sfe,                                         #129-130 
                    slen,                                        #131-132 
                    styp,                                        #133-134 
                    stas,                                        #135-136 
                    stae,                                        #137-138 
                    tatyp,                                       #139-140 
                    afilf,                                       #141-142 
                    afils,                                       #143-144 
                    nofilf,                                      #145-146 
                    nofils,                                      #147-148 
                    lcf,                                         #149-150 
                    hcf,                                         #151-152 
                    lcs,                                         #153-154 
                    hcs,                                         #155-156 
                    year,                                        #157-158 
                    day,                                         #159-160 
                    hour,                                        #161-162 
                    minute,                                      #163-164 
                    sec,                                         #165-166 
                    timbas,                                      #167-168 
                    trwf,                                        #169-170 
                    grnors,                                      #171-172 
                    grnofr,                                      #173-174 
                    grnlof,                                      #175-176 
                    gaps,                                        #177-178 
                    otrav,                                       #179-180 
                    d1,                                          #181-184 
                    f1,                                          #185-188 
                    d2,                                          #189-192 
                    f2,                                          #193-196 
                    ungpow,                                      #197-200
                    unscale,                                     #201-204 
                    ntr,                                         #205-208 
                    mark,                                        #209-210 
                    shortpad,                                    #211-212 
                    unassigned)
        end
end

mutable struct SU_trc
        hdr::SU_hdr
        dat::Array{Float32,1}
        function SU_trc(n::Int64)
                hdr=SU_hdr();
                dat=zeros(Float32,n)
                new(hdr,dat)
        end

end

"""
SUPutTrace(io::IOStream,trc::SU_trc)

Write a SU trace to a file

### Signature
* io  : file stream reference
* trc : seismic trace (header + data)
"""
function SUPutTrace(io::IOStream,trc::SU_trc)
        @assert iswritable(io)
        write(io,Int32(trc.hdr.tracl))
        write(io,Int32(trc.hdr.tracr))
        write(io,Int32(trc.hdr.fldr))
        write(io,Int32(trc.hdr.tracf))
        write(io,Int32(trc.hdr.ep))
        write(io,Int32(trc.hdr.cdp))
        write(io,Int32(trc.hdr.cdpt))
        write(io,Int16(trc.hdr.trid))
        write(io,Int16(trc.hdr.nvs))
        write(io,Int16(trc.hdr.nhs))
        write(io,Int16(trc.hdr.duse))
        write(io,Int32(trc.hdr.offset))
        write(io,Int32(trc.hdr.gelev))
        write(io,Int32(trc.hdr.selev))
        write(io,Int32(trc.hdr.sdepth))
        write(io,Int32(trc.hdr.gdel))
        write(io,Int32(trc.hdr.sdel))
        write(io,Int32(trc.hdr.swdep))
        write(io,Int32(trc.hdr.gwdep))
        write(io,Int16(trc.hdr.scalel))
        write(io,Int16(trc.hdr.scalco))
        write(io,Int32(trc.hdr.sx))
        write(io,Int32(trc.hdr.sy))
        write(io,Int32(trc.hdr.gx))
        write(io,Int32(trc.hdr.gy))
        write(io,Int16(trc.hdr.counit))
        write(io,Int16(trc.hdr.wevel))
        write(io,Int16(trc.hdr.swevel))
        write(io,Int16(trc.hdr.sut))
        write(io,Int16(trc.hdr.gut))
        write(io,Int16(trc.hdr.sstat))
        write(io,Int16(trc.hdr.gstat))
        write(io,Int16(trc.hdr.tstat))
        write(io,Int16(trc.hdr.laga))
        write(io,Int16(trc.hdr.lagb))
        write(io,Int16(trc.hdr.delrt))
        write(io,Int16(trc.hdr.muts))
        write(io,Int16(trc.hdr.mute))
        write(io,Int16(trc.hdr.ns))
        write(io,Int16(trc.hdr.dt))
        write(io,Int16(trc.hdr.gain))
        write(io,Int16(trc.hdr.igc))
        write(io,Int16(trc.hdr.igi))
        write(io,Int16(trc.hdr.corr))
        write(io,Int16(trc.hdr.sfs))
        write(io,Int16(trc.hdr.sfe))
        write(io,Int16(trc.hdr.slen))
        write(io,Int16(trc.hdr.styp))
        write(io,Int16(trc.hdr.stas))
        write(io,Int16(trc.hdr.stae))
        write(io,Int16(trc.hdr.tatyp))
        write(io,Int16(trc.hdr.afilf))
        write(io,Int16(trc.hdr.afils))
        write(io,Int16(trc.hdr.nofilf))
        write(io,Int16(trc.hdr.nofils))
        write(io,Int16(trc.hdr.lcf))
        write(io,Int16(trc.hdr.hcf))
        write(io,Int16(trc.hdr.lcs))
        write(io,Int16(trc.hdr.hcs))
        write(io,Int16(trc.hdr.year))
        write(io,Int16(trc.hdr.day))
        write(io,Int16(trc.hdr.hour))
        write(io,Int16(trc.hdr.minute))
        write(io,Int16(trc.hdr.sec))
        write(io,Int16(trc.hdr.timbas))
        write(io,Int16(trc.hdr.trwf))
        write(io,Int16(trc.hdr.grnors))
        write(io,Int16(trc.hdr.grnofr))
        write(io,Int16(trc.hdr.grnlof))
        write(io,Int16(trc.hdr.gaps))
        write(io,Int16(trc.hdr.otrav))
        write(io,Float32(trc.hdr.d1))
        write(io,Float32(trc.hdr.f1))
        write(io,Float32(trc.hdr.d2))
        write(io,Float32(trc.hdr.f2))
        write(io,Float32(trc.hdr.ungpow))
        write(io,Float32(trc.hdr.unscale))
        write(io,Int32(trc.hdr.ntr))
        write(io,Int16(trc.hdr.mark))
        write(io,Int16(trc.hdr.shortpad))
        @inbounds for i=1:14
                write(io,Int16(trc.hdr.unassigned[i]))
        end
        @inbounds for is=1:trc.hdr.ns
                write(io,Float32(trc.dat[is]))
        end 
end # SUPutTrace

"""
SUPutTrace(io::IOStream,trc::SU_trc)

Read a SU trace from a file

### Signature
* io  : file stream reference
* trc : seismic trace (header + data)

"""
function SUGetTrace(io::IOStream,trc::SU_trc) 
        @assert isreadable(io)
        trc.hdr.tracl=read(io,Int32)
        trc.hdr.tracr=read(io,Int32)
        trc.hdr.fldr=read(io,Int32)
        trc.hdr.tracf=read(io,Int32)
        trc.hdr.ep=read(io,Int32)
        trc.hdr.cdp=read(io,Int32)
        trc.hdr.cdpt=read(io,Int32)
        trc.hdr.trid=read(io,Int16)
        trc.hdr.nvs=read(io,Int16)
        trc.hdr.nhs=read(io,Int16)
        trc.hdr.duse=read(io,Int16)
        trc.hdr.offset=read(io,Int32)
        trc.hdr.gelev=read(io,Int32)
        trc.hdr.selev=read(io,Int32)
        trc.hdr.sdepth=read(io,Int32)
        trc.hdr.gdel=read(io,Int32)
        trc.hdr.sdel=read(io,Int32)
        trc.hdr.swdep=read(io,Int32)
        trc.hdr.gwdep=read(io,Int32)
        trc.hdr.scalel=read(io,Int16)
        trc.hdr.scalco=read(io,Int16)
        trc.hdr.sx=read(io,Int32)
        trc.hdr.sy=read(io,Int32)
        trc.hdr.gx=read(io,Int32)
        trc.hdr.gy=read(io,Int32)
        trc.hdr.counit=read(io,Int16)
        trc.hdr.wevel=read(io,Int16)
        trc.hdr.swevel=read(io,Int16)
        trc.hdr.sut=read(io,Int16)
        trc.hdr.gut=read(io,Int16)
        trc.hdr.sstat=read(io,Int16)
        trc.hdr.gstat=read(io,Int16)
        trc.hdr.tstat=read(io,Int16)
        trc.hdr.laga=read(io,Int16)
        trc.hdr.lagb=read(io,Int16)
        trc.hdr.delrt=read(io,Int16)
        trc.hdr.muts=read(io,Int16)
        trc.hdr.mute=read(io,Int16)
        trc.hdr.ns=read(io,Int16)
        trc.hdr.dt=read(io,Int16)
        trc.hdr.gain=read(io,Int16)
        trc.hdr.igc=read(io,Int16)
        trc.hdr.igi=read(io,Int16)
        trc.hdr.corr=read(io,Int16)
        trc.hdr.sfs=read(io,Int16)
        trc.hdr.sfe=read(io,Int16)
        trc.hdr.slen=read(io,Int16)
        trc.hdr.styp=read(io,Int16)
        trc.hdr.stas=read(io,Int16)
        trc.hdr.stae=read(io,Int16)
        trc.hdr.tatyp=read(io,Int16)
        trc.hdr.afilf=read(io,Int16)
        trc.hdr.afils=read(io,Int16)
        trc.hdr.nofilf=read(io,Int16)
        trc.hdr.nofils=read(io,Int16)
        trc.hdr.lcf=read(io,Int16)
        trc.hdr.hcf=read(io,Int16)
        trc.hdr.lcs=read(io,Int16)
        trc.hdr.hcs=read(io,Int16)
        trc.hdr.year=read(io,Int16)
        trc.hdr.day=read(io,Int16)
        trc.hdr.hour=read(io,Int16)
        trc.hdr.minute=read(io,Int16)
        trc.hdr.sec=read(io,Int16)
        trc.hdr.timbas=read(io,Int16)
        trc.hdr.trwf=read(io,Int16)
        trc.hdr.grnors=read(io,Int16)
        trc.hdr.grnofr=read(io,Int16)
        trc.hdr.grnlof=read(io,Int16)
        trc.hdr.gaps=read(io,Int16)
        trc.hdr.otrav=read(io,Int16)
        trc.hdr.d1=read(io,Float32)
        trc.hdr.f1=read(io,Float32)
        trc.hdr.d2=read(io,Float32)
        trc.hdr.f2=read(io,Float32)
        trc.hdr.ungpow=read(io,Float32)
        trc.hdr.unscale=read(io,Float32)
        trc.hdr.ntr=read(io,Int32)
        trc.hdr.mark=read(io,Int16)
        trc.hdr.shortpad=read(io,Int16)
        trc.hdr.unassigned=zeros(Int16,14)
        for i=1:14
                trc.hdr.unassigned[i]=read(io,Int16)
        end

        if size(trc.dat,1) == undef 
                trc.dat=zeros(trc.hdr.ns)
        end
        @inbounds for is=1:trc.hdr.ns
                trc.dat[is]=read(io,Float32)
        end
        nothing
end # SUGetTrace
end # module SU_mod
