from flask import Flask, request, render_template, Response, jsonify, abort #, make_response, request, render_template, abort, redirect
import os
import sqlite3
from flask.helpers import url_for
import matplotlib
import re
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
def parts_view(part):
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
    transcripts = [ unwrap_transcript([_[8], _[0]] +list(_[9:12])) for _ in data ]
    #print(data)
    #print(data_self)
    #print(tissues)
    #print(transcripts)
    
    return render_template("gene.html", 
                            gene=unwrap_gene(data_self), 
                            transcripts = transcripts, 
                            tissues=tissues)
    #return str(data_self) + "\n######\n"+ str(data)
        
def unwrap_transcript(row):
    """
    Generate a transcript Dictionary from a Transcripts table row
    """
    return {
        'Ensembl_Transcript_ID' : row[0],
        'Ensembl_Gene_ID'       : row[1],
        'Transcript_Start'      : row[2],
        'Transcript_End'        : row[3],
        'Transcript_Biotype'    : row[4]
    }
def unwrap_gene(row):
    """
    Generate a gene Dictionary from a Genes table row
    """
    return {
        'Ensembl_Gene_ID'      : row[0],
        'Associated_Gene_Name' : row[1],
        'Chromosome_Name'      : row[2],
        'Band'                 : row[3],
        'Strand'               : row[4],
        'Gene_Start'           : row[5],
        'Gene_End'             : row[6],
        'Transcript_Count'     : row[7]
    }

@app.route("/transcript/<t_id>")
def transcript_view(t_id):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT * from Transcripts WHERE ensembl_transcript_id = "{t_id}"
    """
    )
    data = [ _ for _ in cur.fetchall() if _ ][0]
    print(data)
    return render_template("transcript.html", transcript = unwrap_transcript(data) )

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
    """ 
    Generate a simple and bugged svg elements
    representing transcripts of a gene
    """
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

       # print(start, gene_start, rel_start)
       # print(end, gene_end, rel_end)
       # print(f"(({rel_end} - {rel_start}) / {gene_len} ) * {_width}")
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
    return rects


### UTILITY VIEWS NOT ASKED FOR IN TP
@app.route("/all_genes")
def view_all_genes():
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute("""
    SELECT * FROM Genes as g
    ORDER BY Transcript_Count
    """)
    data = [ unwrap_gene(_) for _ in cur.fetchall() if _ ]
    #print(data)
    return str(data)

@app.route("/all_transcripts")
def view_all_transcripts():
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute("SELECT * FROM Transcripts GROUP BY ensembl_gene_id")
    data = [ unwrap_transcript(_) for _ in cur.fetchall() if _ ]
    return str(data)




### SECTION TP 7: API-WEB
@app.route('/api/genes/<id>')
def get_gene_detailed(id):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT * from Genes WHERE ensembl_gene_id = "{id}"
    """
    )
    data_self = [ _ for _ in cur.fetchall() if _ ][0]
    if not data_self:
        abort(404)

    detailed_gene = unwrap_gene(data_self)

    cur = con.execute( 
    f"""SELECT *
    FROM Genes as g 
    NATURAL JOIN Transcripts as t
    NATURAL JOIN Expression as e
    WHERE g.ensembl_gene_id = "{id}"
    """
    )

    data = [ _ for _ in cur.fetchall() if _ ]
    detailed_gene["transcripts"] = [ unwrap_transcript([_[8], _[0]] +list(_[9:12])) for _ in data ]
    
    return jsonify(detailed_gene)


def get_gene_compact(id):
    print("->", id)
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT * from Genes WHERE ensembl_gene_id = "{id}"
    """
    )
    data_self = [ _ for _ in cur.fetchall() if _ ][0]
    if not data_self:
        return None
    gene = unwrap_gene(data_self)
    gene["href"] = url_for('get_gene_detailed', id=id)
    return gene

@app.route("/api/genes/", methods=['GET', 'POST'])
def get_many_gene_compact():
    if request.method == 'POST':
        _ = add_gene(request.form)
        return _

    offset = request.args.get('offset', default = 0, type = int)
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
    f"""
    SELECT Ensembl_Gene_ID from Genes limit {offset - 1}, {offset + 100}
    """
    )
    gene_collection = [ get_gene_compact(_[0]) for _ in cur.fetchall() ]
    return jsonify(gene_collection)

def add_gene(post_datum):
    data = []
    print("--->" + str(post_datum))
    
    for k in ['Ensembl_Gene_ID', 'Chromosome_Name', 'Band']:
        if not k in post_datum:
            return jsonify({"error": f"{k} manquant"}), 400
        data.append(post_datum[k])
    for k in ['Gene_Start', 'Gene_End']:
        if k not in post_datum:
            return jsonify({"error": f"{k} manquant"}), 400
        if not type(post_datum[k]) == 'int' :
            if not re.match("[-+]?\d+$", post_datum[k]):
                return jsonify({"error": f"{k} n'est pas un entier"}), 400
        data.append(int(post_datum[k]))
    if not post_datum['Gene_Start'] < post_datum['Gene_End']:
        return jsonify({"error": f"Gene_Start n'est pas inférieure à Gene_End"}), 400
    
    if 'Strand' in post_datum:
        if not type(post_datum['Strand']) == 'int':
            if not re.match("[-+]?\d+$", post_datum['Strand']):
                return jsonify({"error": "Strand n'est pas un entier"}), 400
        data.append(int(post_datum['Strand']))
    else:
        data.append(None)
    if 'Associated_Gene_Name' in post_datum:
        if not type('Associated_Gene_Name') == 'string':
            return jsonify({"error": "Associated_Gene_Name n'est pas une chaine"}), 400
        data.append(post_datum['Associated_Gene_Name'])
    else:
        data.append(None)

    if insert_gene(data):
        return jsonify({ "created": url_for('gene_view',id=data[0]) }), 201
    
    return jsonify({ "error": "Failed to insert" }), 400

def insert_gene(data):
    con = sqlite3.connect("ensembl_hs63_simple.sqlite")
    cur = con.execute( 
        f"""
        INSERT INTO Genes ('Ensembl_Gene_ID', 'Chromosome_Name', 'Band', 'Gene_Start', 'Gene_End', 'Strand', 'Associated_Gene_Name')
        VALUES ('{data[0]}', '{data[1]}', '{data[2]}', '{data[3]}', '{data[4]}', '{data[5]}','{data[6]}')
        """)
    con.commit()
    con.close()
    return True
