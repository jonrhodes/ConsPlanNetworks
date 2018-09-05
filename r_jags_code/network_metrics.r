library('igraph')

#create networks
Star <- graph_from_literal(1-2,1-3,1-4,1-5,1-6)
Wheel <- graph_from_literal(1-2,1-3,1-4,1-5,1-6,2-3,3-4,4-5,5-6,6-2)
Line <- graph_from_literal(1-2,2-3,3-4,4-5,5-6)
Ring <- graph_from_literal(1-2,2-3,3-4,4-5,5-6,6-1)
RingStar <- graph_from_literal(1-2,1-3,2-3,1-4,4-5,4-6)
StarLine <- graph_from_literal(1-2,1-3,1-4,4-5,5-6)
RingLine <- graph_from_literal(1-2,1-3,2-3,1-4,4-5,5-6)
CaseStudy <- graph_from_literal(1-3,1-2,2-3,3-4,5)

#plot case study
plot(CaseStudy,vertex.label=c("Handline","Seine net","Ring net","Gill net","Spear gun"),vertex.size=50,vertex.color="red",vertex.label.font=14,edge.width=3,edge.color="black")

#calculate closeness centrality
CStar <- centralization.closeness(Star)
CWheel <- centralization.closeness(Wheel)
CRing <- centralization.closeness(Ring)
CLine <- centralization.closeness(Line)
CRingStar <- centralization.closeness(RingStar)
CStarLine <- centralization.closeness(StarLine)
CRingLine <- centralization.closeness(RingLine)
CCaseStudy <- centralization.closeness(CaseStudy)

CStar
CWheel
CRing
CLine
CRingStar
CStarLine
CRingLine
CCaseStudy

#calculate degree centrality
DStar <- centralization.degree(Star)
DWheel <- centralization.degree(Wheel)
DRing <- centralization.degree(Ring)
DLine <- centralization.degree(Line)
DRingStar <- centralization.degree(RingStar)
DStarLine <- centralization.degree(StarLine)
DRingLine <- centralization.degree(RingLine)
DCaseStudy <- centralization.degree(CaseStudy)

DStar
DWheel
DRing
DLine
DRingStar
DStarLine
DRingLine
DCaseStudy

#calculate betweenness centrality
BStar <- centralization.betweenness(Star)
BWheel <- centralization.betweenness(Wheel)
BRing <- centralization.betweenness(Ring)
BLine <- centralization.betweenness(Line)
BRingStar <- centralization.betweenness(RingStar)
BStarLine <- centralization.betweenness(StarLine)
BRingLine <- centralization.betweenness(RingLine)
BCaseStudy <- centralization.betweenness(CaseStudy)

BStar
BWheel
BRing
BLine
BRingStar
BStarLine
BRingLine
BCaseStudy

#calculare graph density
DensStar <- graph.density(Star)
DensWheel <- graph.density(Wheel)
DensRing <- graph.density(Ring)
DensLine <- graph.density(Line)
DensRingStar <- graph.density(RingStar)
DensStarLine <- graph.density(StarLine)
DensRingLine <- graph.density(RingLine)
DensCaseStudy <- graph.density(CaseStudy)

DensStar
DensWheel
DensRing
DensLine
DensRingStar
DensStarLine
DensRingLine
DensCaseStudy





