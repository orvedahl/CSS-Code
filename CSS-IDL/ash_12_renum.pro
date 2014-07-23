;+
;     ash_12_renum
;     Ben Brown
;     
;     Renumbers files to match ASH_1.2 7-digit numbering.  This is
;     done by appending padding zeros to the beginning of numbers
;     smaller than 1.0d6
;-
function ash_12_renum, files, ash_12_numbering_starts=ash_12_numbering_starts, quiet=quiet

    ash_digits = 7
    if not(keyword_set(ash_12_numbering_starts)) then ash_12_numbering_starts = 1

    textnum=strcompress(files,/r)
                 
    text_index = where((files ge ash_12_numbering_starts) and (files lt 1000000.), $
                       text_count) 
    if text_count ne 0 then begin
        for i=0, n_elements(text_index)-1 do begin
            while (strlen(textnum[text_index[i]]) lt ash_digits) do $
                     textnum(text_index[i]) = '0'+textnum(text_index[i])
        endfor
        if not keyword_set(quiet) then print, 'converting: ', textnum(text_index)
    endif
    return, textnum
end
