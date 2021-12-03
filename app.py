from flask import Flask, render_template, Response #, make_response, request, render_template, abort, redirect
import os
import sqlite3
import matplotlib

app = Flask(__name__)

@app.route('/plot.svg')
def get_svg():
    return render_template("generate_svg.html")
@app.route('/plot2.svg')
def get_svg2():
    return Response(
        """
        <svg
            xmlns="http://www.w3.org/2000/svg" 
            xml:lang="en" 
            xmlns:xlink="http://www.w3.org/1999/xlink"
            viewBox="-50 -50 100 100" >
                <circle cx="0" cy="0" r="40" style="fill:red" />
                <rect x="-27.5" y="-7.5" width="55" height="15"
                    style="fill:white" />
        </svg>
    """,
        mimetype='image/svg+xml'
    )
@app.route("/")
def home():
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    """
    SELECT DISTINCT e.atlas_organism_part
    FROM Expression as e
    ORDER BY atlas_organism_part ASC
    """
    )
    data = [ _[0] for _ in cur.fetchall() if _[0] ]
    return render_template("home.html", tissues=data)

@app.route("/parts/<part>/genes")
def view_part(part):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT DISTINCT g.ensembl_gene_id, associated_gene_name
    FROM Genes as g
    NATURAL JOIN Transcripts as t
    NATURAL JOIN Expression as e
    WHERE atlas_organism_part = "{part}"
    ORDER BY g.ensembl_gene_id
    """
    )
    data = [ _ for _ in cur.fetchall() if _ ]
  
    return render_template("parts.html", data=data)

@app.route("/gene/<id>")
def gene_view(id):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT * from Genes WHERE ensembl_gene_id = "{id}"
    """
    )
    data_self = [ _ for _ in cur.fetchall() if _ ][0]
    cur = con.execute( 
    f"""SELECT *
    FROM Genes as g 
    NATURAL JOIN Transcripts as t
    NATURAL JOIN Expression as e
    WHERE g.ensembl_gene_id = "{id}"
    """
    )
    data = [ _ for _ in cur.fetchall() if _ ]
    tissues = list( set([_[-1] for _ in data if _[-1]]) )
    transcripts = [ _[8:12] for _ in data ]

    print(data_self)
    print(tissues)
    print(transcripts)
    
    return render_template("gene.html", 
                            gene=data_self, 
                            transcripts = transcripts, 
                            tissues=tissues)
    #return str(data_self) + "\n######\n"+ str(data)

class Transcript():
    def __init__(self, *args):
        self.Ensembl_Transcript_ID = args[0]
        self.Ensembl_Gene_ID       = args[1]
        self.Transcript_Start      = args[2]
        self.Transcript_End        = args[3]
        self.Transcript_Biotype    = args[4]

@app.route("/transcripts/<gene>/histogram.svg")
def get_transcripts_per_gene_svg(gene):
   
    w = 500
    h = 500

    svg_string = f"""
        <svg
            xmlns="http://www.w3.org/2000/svg" 
            xml:lang="en" 
            xmlns:xlink="http://www.w3.org/1999/xlink"
            viewBox="0 0 {w} {h}" >
    """
    rects =  get_transcripts_svg_data(gene, w, h)
    for rec_x, rec_y, rec_w, rec_h, t_id, hex in rects:
       # rgb_str = ','.join([str(_) for _ in rgba])
        svg_string += f"""<rect x="{rec_x}" y="{rec_y}" 
                        width="{rec_w}" height="{rec_h}" 
                        style='fill:{hex}'/>\n"""

    svg_string += "</svg>"
    return Response(
        svg_string
        , mimetype='image/svg+xml'
    )

def get_transcripts_svg_data(gene, width  = 500, height = 500):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    """
    SELECT DISTINCT g.ensembl_gene_id, g.associated_gene_name, g.gene_start, g.gene_end, t.ensembl_transcript_id, t.transcript_start, t.transcript_end, t.transcript_biotype
    FROM Genes as g 
    NATURAL JOIN Transcripts as t
    WHERE g.ensembl_gene_id = '""" + gene + """'
    """
    )
    data = [ _ for _ in cur.fetchall() if _ ]
    
    margin = 10
    
    cst_rec_height = (height - 2 * margin) / len(data)

    
    def xscale(start, end):
        gene_start = data[0][2]
        gene_end  = data[0][3] 
        gene_len = gene_end - gene_start

        rel_start = start - gene_start
        rel_end   = end - gene_start
        _width = width - margin * 2

        print(start, gene_start, rel_start)
        print(end, gene_end, rel_end)
        print(f"(({rel_end} - {rel_start}) / {gene_len} ) * {_width}")
        svg_width = ((rel_end - rel_start) / gene_len ) * _width
        svg_start = margin + ( rel_start / gene_len ) * _width
        return svg_start, svg_width

    curr_y = margin
    rects = []
    cmap = matplotlib.cm.get_cmap('tab20c')

    for i, (g_id, g_name, g_start, g_end, t_id, t_start, t_end, t_biotype) in enumerate(data):
        x_start, rec_len  = xscale(t_start, t_end)
        rects.append(
           (x_start, margin + i * cst_rec_height,
            rec_len, cst_rec_height,
            t_id, matplotlib.colors.to_hex( cmap(i) ))
        )
    print(rects)
    return rects

@app.route("/all_genes")
def root():
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute("""
    SELECT g.ensembl_gene_id, g.associated_gene_name, COUNT(*) AS transcript_count 
    FROM Genes as g
    LEFT JOIN Transcripts as t ON t.ensembl_gene_id = g.ensembl_gene_id
    GROUP BY 1
    """)
    data = [ _ for _ in cur.fetchall() if _ ]
    #print(data)
    return str(data)



@app.route("/all_transcript")
def transcripts():
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute("SELECT * FROM Transcripts GROUP BY ensembl_gene_id")
    data = [ _ for _ in cur.fetchall() if _ ]
    return str(data)