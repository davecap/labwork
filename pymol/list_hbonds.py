from pymol import cmd

def list_hb(selection,cutoff=3.2,angle=55,hb_list_name='hbonds'):
  """
  USAGE

  list_hb selection, [cutoff (default=3.2)], [angle (default=55)], [hb_list_name]
  
  e.g.
    list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds
  """
  cutoff=float(cutoff)
  angle=float(angle)
  hb = cmd.find_pairs("((byres "+selection+") and n;n)","((byres "+selection+") and n;o)",mode=1,cutoff=cutoff,angle=angle)
# sort the list for easier reading
  hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))

  for pairs in hb:
    cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print "%s/%3s`%s/%s " % (chain,resn,resi,name),')
    cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print "%s/%3s`%s/%s " % (chain,resn,resi,name),')
    print "%.2f" % cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))

cmd.extend("list_hb",list_hb)
