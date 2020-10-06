def bold(t): 
    t = f"""<b>{t}</b>""" 
    return(t)

def italic(t):
    t = f"""<i>{t}</i>"""
    return(t)

def under_line(t):
    t = f"""<u>{t}</u>"""
    return(t)

def tag_color(t,c):
    t = f"""<font color={c}>{t}</font>"""
    return(t)

def purple(t):
    t = tag_color(t, "#0000ff")
    return(t)

def red(t):
    t = tag_color(t, "#AA0000")
    return(t)


def fontsize(t,size="2.5"):
    t = f"""<font size={size}>{t}</font>"""
    #t = f"""<span style="font-size:{fontsize}"> {t} </span>"""
    return(t)
