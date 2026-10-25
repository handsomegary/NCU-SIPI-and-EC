

def snpreader(filename,maxsamplenum=5000,reciproc=True):
    ''' Read touchstone .snp files 
    Input: .snp file, maxsamplenum
        maxsamplenum: maximum # of freq samples
        reciproc: reciprocal network. Read only the 
            upper triangular part of the matrix if True.
    Output: freq, sp, Z0
        freq: frequency vector. Nfreq x 1
        sp: S parameters in numpy array. Nfreq x Nc
    
    Example:
    
    from snpreader import snpreader
    fname = 'your_file.s2p'
    reci = True
    maxnum = 5000
    data = snpreader(fname,reciproc=reci,maxsamplenum=maxnum)
    freq = data[0]
    sp = data[1]
    Z0 = data[2]
    
    # sp column indexing example (4-port, reciproc=False)
    # Column Sij
    #    0   S11
    #    1   S12
    #    2   S13
    #    3   S14
    #    4   S21
    #    5   S22
    #    6   S23
    #    7   S24
    #    8   S31
    #    9   S32
    #   10   S33
    #   11   S34
    #   12   S41
    #   13   S42
    #   14   S43
    #   15   S44
    
    # sp column indexing example (4-port, reciproc=True)
    # Column Sij
    #    0   S11
    #    1   S12
    #    2   S13
    #    3   S14
    #    4   S22
    #    5   S23
    #    6   S24
    #    7   S33
    #    8   S34
    #    9   S44
    
    
    
    The reading is in a 2-step manner. The first
    run reads the frequency and obtain the total
    number of freq points. In the second run, the
    data points are extracted, with possibly a
    reduced sample number according to the input
    'maxsamplenum'.
    Please renormalize all ports and all 
    frequencies to the same real inpedance,
    as specified in the option line.'''
    
    Nport = 0
    for i in range(10000):
        if '.s'+str(i+1)+'p' in filename.lower():
            Nport = i+1
            break
    if Nport==0:
        print('Please include .snp in the filename,')
        print('where n is the number of ports.')
        return None
    print('Nport = %d'%Nport)
    print('Read frequency data... ', end='\n')
    
    import numpy as np
    from time import time
    tic = time()
    Ntemp=3999 # could be any large number
    freq=np.zeros(1)
    
    # First run: identify option line and number of points
    # with open(filename,mode='r',encoding='utf-8') as f:
    with open(filename,mode='r') as f:
        for line in f:
            if line[0]=='#':
                optionline=line
                break
        idx=optionline.find('!') # in-line comment
        if idx!=-1:
            optionline=optionline[:idx] # truncate comments
        idxf=0
        count=0
        numele = 2*(Nport**2)
        for line in f:
            if line[0]=='!': # ignore comment line
                continue
            if line[0]=='#': # ignore multiple option line
                continue
            idx=line.find('!') # in-line comment
            if idx!=-1:
                line=line[:idx] # truncate comments
            content=line.split()
            if not content: # it is an empty line
                continue
            #!!!
            if count == 0:
                if idxf % Ntemp == 0:
                    freq = np.hstack((freq,
                            np.zeros(Ntemp)))
                freq[idxf]=float(content[0])
            count += len(content)
            if count > numele:
                count = 0
                idxf += 1
            
            # for item in content:
            #     if count == 0:
            #         if idxf % Ntemp == 0:
            #             freq = np.hstack((freq,
            #                     np.zeros(Ntemp)))
            #         count += 1
            #         freq[idxf]=float(item)
            #         continue
            #     count += 1
            #     if count > numele:
            #         count = 0
            #         idxf += 1
    Nfreq0 = idxf
    freq = freq[:Nfreq0]
    
    # Reduce sample
    lowfs = 0 # number of consecutive low freq points kept
    if (Nfreq0 <= maxsamplenum):
        Nfreq = Nfreq0
        freqidx = np.arange(Nfreq0)
        print('Read all %d frequency points.'%Nfreq,end='\n')
    else:
        print('Read only %d out of %d freq points.'
              %(maxsamplenum,Nfreq0),end='\n')
        Nfreq = maxsamplenum
        # note: the way of reducing samples deserves study
        # Where more points are sampled, where less error can
        # be expected after least-squares fitting
        freqidx_low = np.arange(lowfs,dtype=int)
        freqidx_high = np.around(np.linspace
                    (lowfs,Nfreq0-1,Nfreq-lowfs)).astype(int)
        freqidx = np.hstack((freqidx_low,freqidx_high))
        freq = freq[freqidx]
    if reciproc:
        Nc = (Nport*(Nport+1))//2
        print('Read only the upper triangular part.')
    else:
        Nc = Nport**2
        print('Read the full matrix.')
    # data format
    if 'db' in optionline.lower():
        dataformat = 0
        print('db ', end='')
    elif 'ma' in optionline.lower():
        dataformat = 1
        print('ma ',end='')
    elif 'ri' in optionline.lower():
        dataformat = 2
        print('ri ',end='')
    else:
        dataformat = 1
        print('ma (default) ')
        
    
    # frequency unit
    if 'ghz' in optionline.lower():
        multiplier = 10**9
        print('GHz ',end='')
    elif 'mhz' in optionline.lower():
        multiplier = 10**6
        print('MHz ',end='')
    elif 'khz' in optionline.lower():
        multiplier = 10**3
        print('kHz ',end='')
    elif 'hz' in optionline.lower():
        multiplier = 1
        print('Hz ',end='')
    else:
        multiplier = 1e9
        print('GHz (default) ',end='')
    freq = freq * multiplier
    
    # reference impedance
    # Z0 = optionline.split()[-1]
    # try:
    #     Z0 = float(Z0)
    # except ValueError:
    #     Z0 = 50
    try:
        options = optionline.lower().split()
        Z0idx = options.index('r')
        Z0 = float(options[Z0idx+1])
        print('Z0 = %.2f Ohm' % Z0)
    except:
        Z0 = 50
        print('Z0 = 50 Ohm (default)')
    
    # The following commented version is simpler, but may 
    # cost more memory. Sometimes, for >100-port snp files,
    # Python is able to allocate memory for the variables
    # x and y, but fails to allocate memory for sp in the end.
    # Hence, in the new version, I only define sp but not xy.
    # x=np.zeros((Nfreq,Nport**2))
    # y=np.zeros((Nfreq,Nport**2))
    # with open(filename,mode='r',encoding='utf-8') as f:
    #     for line in f:
    #         if line[0]=='#':
    #             break
    #     idxf0 = idxf = 0
    #     nextf = freqidx[idxf]
    #     count = 0
    #     for line in f:
    #         if line[0]=='!': # ignore comment line
    #             continue
    #         idx=line.find('!') # in-line comment
    #         if idx!=-1:
    #             line=line[:idx] # truncate comments
    #         content=line.split()
    #         if not content: # it is an empty line
    #             continue
    #         for item in content:
    #             # if idxf0 not in freqidx:
    #             if idxf0 != nextf:
    #                 count += 1
    #                 if count > numele:
    #                     count = 0
    #                     idxf0 += 1
    #                 continue
    #             if count == 0:
    #                 count += 1
    #                 continue
    #             if count%2 == 1:
    #                 x[idxf,(count//2)]=float(item)
    #                 count += 1
    #             elif count%2 == 0:
    #                 y[idxf,count//2-1]=float(item)
    #                 count += 1
    #             if count > numele:
    #                 count = 0
    #                 idxf += 1
    #                 idxf0 += 1
    #                 if idxf < Nfreq:
    #                     nextf = freqidx[idxf]
    #     # data format
    #     if dataformat == 0:
    #         sp = (10**(x/20)) * np.exp(1j*y/180*np.pi)
    #     elif dataformat == 1:
    #         sp = x * np.exp(1j*y/180*np.pi)
    #     elif dataformat == 2:
    #         sp = x + 1j*y
    
    # New version, more memory efficient
    Nfreq10=Nfreq//10
    print('Read S parameter data: ',end='')
    sp = np.zeros((Nfreq,Nc),dtype=complex)
    spx = np.zeros((Nc,),dtype=float) #!!!
    spy = np.zeros((Nc,),dtype=float) #!!!
    if reciproc:
        # with open(filename,mode='r',encoding='utf-8') as f:
        with open(filename,mode='r') as f:
            for line in f:
                if line[0]=='#':
                    break
            idxf0 = idxf = 0
            nextf = freqidx[idxf]
            i = j = nc = count = 0
            for line in f:
                if line[0]=='!': # ignore comment line
                    continue
                if line[0]=='#': # ignore multiple option line
                    continue
                idx=line.find('!') # in-line comment
                if idx!=-1:
                    line=line[:idx] # truncate comments
                content=line.split()
                if not content: # it is an empty line
                    continue
                for item in content:
                    # if idxf0 not in freqidx:
                    if idxf0 != nextf:
                        count += 1
                        if count > numele:
                            count = 0
                            idxf0 += 1
                        continue
                    if count == 0:
                        count += 1
                        continue
                    if i > j:
                        if count%2 == 0: j += 1
                        count += 1
                        continue
                    if count%2 == 1:
                        # x = float(item)
                        spx[nc] = float(item) #!!!
                        count += 1
                    elif count%2 == 0:
                        # y = float(item)
                        spy[nc] = float(item) #!!!
                        # if dataformat == 0:
                        #     sp[idxf,nc] = (10**(x/20))\
                        #          * np.exp(1j*y/180*np.pi)
                        # elif dataformat == 1:
                        #     sp[idxf,nc] = \
                        #         x * np.exp(1j*y/180*np.pi)
                        # elif dataformat == 2:
                        #     sp[idxf,nc] = x + 1j*y
                        count += 1
                        nc += 1
                        j += 1
                        if j%Nport == 0:
                            i += 1
                            j = 0
                    if count > numele:
                        #!!!
                        if dataformat == 0:
                            sp[idxf,:] = (10**(spx/20)) * np.exp(1j*spy/180*np.pi)
                        elif dataformat == 1:
                            sp[idxf,:] = spx * np.exp(1j*spy/180*np.pi)
                        elif dataformat == 2:
                            sp[idxf,:] = spx + 1j*spy
                        
                        count = nc = i = j = 0
                        idxf += 1
                        idxf0 += 1
                        if idxf < Nfreq:
                            nextf = freqidx[idxf]
                        if idxf%(Nfreq10) == 0:
                            print('%d'%(idxf//(Nfreq10)),
                                                  end='')
    else:
        # with open(filename,mode='r',encoding='utf-8') as f:
        with open(filename,mode='r') as f:
            for line in f:
                if line[0]=='#':
                    break
            idxf0 = idxf = 0
            nextf = freqidx[idxf]
            count = 0
            for line in f:
                if line[0]=='!': # ignore comment line
                    continue
                if line[0]=='#': # ignore multiple option line
                    continue
                idx=line.find('!') # in-line comment
                if idx!=-1:
                    line=line[:idx] # truncate comments
                content=line.split()
                if not content: # it is an empty line
                    continue
                for item in content:
                    # if idxf0 not in freqidx:
                    if idxf0 != nextf:
                        count += 1
                        if count > numele:
                            count = 0
                            idxf0 += 1
                        continue
                    if count == 0:
                        count += 1
                        continue
                    if count%2 == 1:
                        # x = float(item)
                        spx[count//2] = float(item) #!!!
                        count += 1
                    elif count%2 == 0:
                        # y = float(item)
                        spy[count//2-1] = float(item) #!!!
                        # if dataformat == 0:
                        #     sp[idxf,count//2-1] = (10**(x/20))\
                        #          * np.exp(1j*y/180*np.pi)
                        # elif dataformat == 1:
                        #     sp[idxf,count//2-1] = \
                        #         x * np.exp(1j*y/180*np.pi)
                        # elif dataformat == 2:
                        #     sp[idxf,count//2-1] = x + 1j*y
                        count += 1
                    if count > numele:
                        #!!!
                        if dataformat == 0:
                            sp[idxf,:] = (10**(spx/20)) * np.exp(1j*spy/180*np.pi)
                        elif dataformat == 1:
                            sp[idxf,:] = spx * np.exp(1j*spy/180*np.pi)
                        elif dataformat == 2:
                            sp[idxf,:] = spx + 1j*spy
                        
                        count = 0
                        idxf += 1
                        idxf0 += 1
                        if idxf < Nfreq:
                            nextf = freqidx[idxf]
                        if idxf%(Nfreq10) == 0:
                            print('%d'%(idxf//(Nfreq10)),
                                                  end='')
        if Nport == 2:
            temp = sp[:,1].copy()
            sp[:,1] = sp[:,2]
            sp[:,2] = temp
    
    print('\nCompleted in %d sec.'%(time()-tic))
    return freq, sp, Z0




# def snpreader(filename,maxsamplenum=5000,reciproc=True):
#     ''' Read touchstone .snp files 
#     Input: .snp file, maxsamplenum
#         maxsamplenum: maximum # of freq samples
#         reciproc: reciprocal network. Read only the 
#             upper triangular part of the matrix if True.
#     Output: freq, sp, Z0
#         freq: frequency vector. Nfreq x 1
#         sp: S parameters in numpy array. Nfreq x Nc
    
#     Example:
    
#     from snpreader import snpreader
#     fname = 'your_file.s2p'
#     reci = True
#     maxnum = 5000
#     data = snpreader(fname,reciproc=reci,maxsamplenum=maxnum)
#     freq = data[0]
#     sp = data[1]
#     Z0 = data[2]
    
    
#     The reading is in a 2-step manner. The first
#     run reads the frequency and obtain the total
#     number of freq points. In the second run, the
#     data points are extracted, with possibly a
#     reduced sample number according to the input
#     'maxsamplenum'.
#     Please renormalize all ports and all 
#     frequencies to the same real inpedance,
#     as specified in the option line.'''
    
#     Nport = 0
#     for i in range(10000):
#         if '.s'+str(i+1)+'p' in filename.lower():
#             Nport = i+1
#             break
#     if Nport==0:
#         print('Please include .snp in the filename,')
#         print('where n is the number of ports.')
#         return None
#     print('Nport = %d'%Nport)
#     print('Read frequency data... ', end='\n')
    
#     import numpy as np
#     from time import time
#     tic = time()
#     Ntemp=3999 # could be any large number
#     freq=np.zeros(1)
    
#     # First run: identify option line and number of points
#     with open(filename,mode='r',encoding='utf-8') as f:
#         for line in f:
#             if line[0]=='#':
#                 optionline=line
#                 break
#         idx=optionline.find('!') # in-line comment
#         if idx!=-1:
#             optionline=optionline[:idx] # truncate comments
#         idxf=0
#         count=0
#         numele = 2*(Nport**2)
#         for line in f:
#             if line[0]=='!': # ignore comment line
#                 continue
#             if line[0]=='#': # ignore multiple option line
#                 continue
#             idx=line.find('!') # in-line comment
#             if idx!=-1:
#                 line=line[:idx] # truncate comments
#             content=line.split()
#             if not content: # it is an empty line
#                 continue
#             for item in content:
#                 if count == 0:
#                     if idxf % Ntemp == 0:
#                         freq = np.hstack((freq,
#                                 np.zeros(Ntemp)))
#                     count += 1
#                     freq[idxf]=float(item)
#                     continue
#                 count += 1
#                 if count > numele:
#                     count = 0
#                     idxf += 1
#     Nfreq0 = idxf
#     freq = freq[:Nfreq0]
    
#     # Reduce sample
#     lowfs = 0 # number of consecutive low freq points kept
#     if (Nfreq0 <= maxsamplenum):
#         Nfreq = Nfreq0
#         freqidx = np.arange(Nfreq0)
#         print('Read all %d frequency points.'%Nfreq,end='\n')
#     else:
#         print('Read only %d out of %d freq points.'
#               %(maxsamplenum,Nfreq0),end='\n')
#         Nfreq = maxsamplenum
#         # note: the way of reducing samples deserves study
#         # Where more points are sampled, where less error can
#         # be expected after least-squares fitting
#         freqidx_low = np.arange(lowfs,dtype=int)
#         freqidx_high = np.around(np.linspace
#                     (lowfs,Nfreq0-1,Nfreq-lowfs)).astype(int)
#         freqidx = np.hstack((freqidx_low,freqidx_high))
#         freq = freq[freqidx]
#     if reciproc:
#         Nc = (Nport*(Nport+1))//2
#         print('Read only the upper triangular part.')
#     else:
#         Nc = Nport**2
#         print('Read the full matrix.')
#     # data format
#     if 'db' in optionline.lower():
#         dataformat = 0
#         print('db ', end='')
#     elif 'ma' in optionline.lower():
#         dataformat = 1
#         print('ma ',end='')
#     elif 'ri' in optionline.lower():
#         dataformat = 2
#         print('ri ',end='')
#     else:
#         dataformat = 1
#         print('ma (default)',end='')
    
#     # frequency unit
#     if 'ghz' in optionline.lower():
#         multiplier = 10**9
#         print('GHz ',end='')
#     elif 'mhz' in optionline.lower():
#         multiplier = 10**6
#         print('MHz ',end='')
#     elif 'khz' in optionline.lower():
#         multiplier = 10**3
#         print('KHz ',end='')
#     elif 'hz' in optionline.lower():
#         multiplier = 1
#         print('Hz ',end='')
#     else:
#         multiplier = 10**9
#         print('GHz (default)',end='')
#     freq = freq * multiplier
    
#     # reference impedance
#     # Z0 = optionline.split()[-1]
#     # try:
#     #     Z0 = float(Z0)
#     # except ValueError:
#     #     Z0 = 50
#     try:
#         options = optionline.lower().split()
#         Z0idx = options.index('r')
#         Z0 = float(options[Z0idx+1])
#         print('Z0 = %.2f Ohm' % Z0)
#     except:
#         Z0 = 50
#         print('Z0 = 50 Ohm (default)')
    
#     # The following commented version is simpler, but may 
#     # cost more memory. Sometimes, for >100-port snp files,
#     # Python is able to allocate memory for the variables
#     # x and y, but fails to allocate memory for sp in the end.
#     # Hence, in the new version, I only define sp but not xy.
#     # x=np.zeros((Nfreq,Nport**2))
#     # y=np.zeros((Nfreq,Nport**2))
#     # with open(filename,mode='r',encoding='utf-8') as f:
#     #     for line in f:
#     #         if line[0]=='#':
#     #             break
#     #     idxf0 = idxf = 0
#     #     nextf = freqidx[idxf]
#     #     count = 0
#     #     for line in f:
#     #         if line[0]=='!': # ignore comment line
#     #             continue
#     #         idx=line.find('!') # in-line comment
#     #         if idx!=-1:
#     #             line=line[:idx] # truncate comments
#     #         content=line.split()
#     #         if not content: # it is an empty line
#     #             continue
#     #         for item in content:
#     #             # if idxf0 not in freqidx:
#     #             if idxf0 != nextf:
#     #                 count += 1
#     #                 if count > numele:
#     #                     count = 0
#     #                     idxf0 += 1
#     #                 continue
#     #             if count == 0:
#     #                 count += 1
#     #                 continue
#     #             if count%2 == 1:
#     #                 x[idxf,(count//2)]=float(item)
#     #                 count += 1
#     #             elif count%2 == 0:
#     #                 y[idxf,count//2-1]=float(item)
#     #                 count += 1
#     #             if count > numele:
#     #                 count = 0
#     #                 idxf += 1
#     #                 idxf0 += 1
#     #                 if idxf < Nfreq:
#     #                     nextf = freqidx[idxf]
#     #     # data format
#     #     if dataformat == 0:
#     #         sp = (10**(x/20)) * np.exp(1j*y/180*np.pi)
#     #     elif dataformat == 1:
#     #         sp = x * np.exp(1j*y/180*np.pi)
#     #     elif dataformat == 2:
#     #         sp = x + 1j*y
    
#     # New version, more memory efficient
#     Nfreq10=Nfreq//10
#     print('Read S parameter data: ',end='')
#     sp = np.zeros((Nfreq,Nc),dtype=complex)
#     if reciproc:
#         with open(filename,mode='r',encoding='utf-8') as f:
#             for line in f:
#                 if line[0]=='#':
#                     break
#             idxf0 = idxf = 0
#             nextf = freqidx[idxf]
#             i = j = nc = count = 0
#             for line in f:
#                 if line[0]=='!': # ignore comment line
#                     continue
#                 if line[0]=='#': # ignore multiple option line
#                     continue
#                 idx=line.find('!') # in-line comment
#                 if idx!=-1:
#                     line=line[:idx] # truncate comments
#                 content=line.split()
#                 if not content: # it is an empty line
#                     continue
#                 for item in content:
#                     # if idxf0 not in freqidx:
#                     if idxf0 != nextf:
#                         count += 1
#                         if count > numele:
#                             count = 0
#                             idxf0 += 1
#                         continue
#                     if count == 0:
#                         count += 1
#                         continue
#                     if i > j:
#                         if count%2 == 0: j += 1
#                         count += 1
#                         continue
#                     if count%2 == 1:
#                         x = float(item)
#                         count += 1
#                     elif count%2 == 0:
#                         y = float(item)
#                         if dataformat == 0:
#                             sp[idxf,nc] = (10**(x/20))\
#                                  * np.exp(1j*y/180*np.pi)
#                         elif dataformat == 1:
#                             sp[idxf,nc] = \
#                                 x * np.exp(1j*y/180*np.pi)
#                         elif dataformat == 2:
#                             sp[idxf,nc] = x + 1j*y
#                         count += 1
#                         nc += 1
#                         j += 1
#                         if j%Nport == 0:
#                             i += 1
#                             j = 0
#                     if count > numele:
#                         count = nc = i = j = 0
#                         idxf += 1
#                         idxf0 += 1
#                         if idxf < Nfreq:
#                             nextf = freqidx[idxf]
#                         if idxf%(Nfreq10) == 0:
#                             print('%d'%(idxf//(Nfreq10)),
#                                                   end='')
#     else:
#         with open(filename,mode='r',encoding='utf-8') as f:
#             for line in f:
#                 if line[0]=='#':
#                     break
#             idxf0 = idxf = 0
#             nextf = freqidx[idxf]
#             count = 0
#             for line in f:
#                 if line[0]=='!': # ignore comment line
#                     continue
#                 if line[0]=='#': # ignore multiple option line
#                     continue
#                 idx=line.find('!') # in-line comment
#                 if idx!=-1:
#                     line=line[:idx] # truncate comments
#                 content=line.split()
#                 if not content: # it is an empty line
#                     continue
#                 for item in content:
#                     # if idxf0 not in freqidx:
#                     if idxf0 != nextf:
#                         count += 1
#                         if count > numele:
#                             count = 0
#                             idxf0 += 1
#                         continue
#                     if count == 0:
#                         count += 1
#                         continue
#                     if count%2 == 1:
#                         x = float(item)
#                         count += 1
#                     elif count%2 == 0:
#                         y = float(item)
#                         if dataformat == 0:
#                             sp[idxf,count//2-1] = (10**(x/20))\
#                                  * np.exp(1j*y/180*np.pi)
#                         elif dataformat == 1:
#                             sp[idxf,count//2-1] = \
#                                 x * np.exp(1j*y/180*np.pi)
#                         elif dataformat == 2:
#                             sp[idxf,count//2-1] = x + 1j*y
#                         count += 1
#                     if count > numele:
#                         count = 0
#                         idxf += 1
#                         idxf0 += 1
#                         if idxf < Nfreq:
#                             nextf = freqidx[idxf]
#                         if idxf%(Nfreq10) == 0:
#                             print('%d'%(idxf//(Nfreq10)),
#                                                   end='')
#         if Nport == 2:
#             temp = sp[:,1].copy()
#             sp[:,1] = sp[:,2]
#             sp[:,2] = temp 
 
#     print('\nCompleted in %d sec.\n'%(time()-tic))
#     return freq, sp, Z0

