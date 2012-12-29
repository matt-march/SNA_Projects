# see http://networkx.lanl.gov/ for a tutorial on NetworkX
import networkx as nx
import numpy as np
import types
from datetime import date
from datetime import datetime
from copy import deepcopy
from pylab import *
from collections import deque

#######################################
#
# Helper functions used
#
#######################################

def makeDate(x): 
    xSplit = x.split('-')
    return date(int(xSplit[0]), int(xSplit[1]), int(xSplit[2]))

def noAttribute(x): return x[1] == {}

def getKeys(x): return x[0]

#######################################
#
# Create the graph
#
#######################################

# Create a graph
G=nx.DiGraph()
G.nodes(data=True)

print 'Building initial graph'

# Open the file
citationsFile = open('C:\Users\matt\SkyDrive\Social Network Analysis - Coursera\Project/Cit-HepPh.txt', 'r')

# Add the edges
citationsLine = citationsFile.readline()
while citationsLine != '':
    splitCitationsLine = citationsLine.split('\t')
    
    try: 
        G.add_edge(int(splitCitationsLine[0]), int(splitCitationsLine[1].rstrip()))
    except ValueError as e: 
        # Do nothing, just means we have a comment, not a real entry
        citationsLine = citationsFile.readline()
        
    citationsLine = citationsFile.readline()

citationsFile.close()

print 'There are ' + str(len(G.nodes())) + ' nodes initially'

# Remove self-loops
# There were some self-loops in the network, so remove them
G.remove_edges_from(G.selfloop_edges())

# Add the dates, track the oldest and newest
print 'Adding dates'
datesFile = open('C:\Users\matt\SkyDrive\Social Network Analysis - Coursera\Project/Cit-HepPh-dates.txt', 'r')

oldestDate = date.max
newestDate = date.min

datesLine = datesFile.readline()
while datesLine != '':
    splitDates = datesLine.split('\t')
    try: 
        dateNode = int(splitDates[0])
        publishedDateString = splitDates[1].rstrip()
        publishedDate = makeDate(publishedDateString)
        G.node[dateNode]['published'] = publishedDateString
        if (publishedDate < oldestDate):
            oldestDate = publishedDate
        if (publishedDate > newestDate):
            newestDate = publishedDate
        datesLine = datesFile.readline()
    except KeyError as e:
        # Do nothing, just means we have a date with no node in the graph
        datesLine = datesFile.readline()
    except ValueError as e:
        # Do nothing, just means we have a comment, not a real entry
        datesLine = datesFile.readline()

datesFile.close()

# Filter out the nodes with no publication date and remove from the graph
removeList = map(getKeys, filter(noAttribute, G.nodes(data=True)))

print 'Removing ' + str(len(removeList)) + ' no-date nodes'
print
G.remove_nodes_from(removeList)

#######################################
#
# Properties of the final graph
#
#######################################

print
print 'Final graph properties'
print

# Node and edge counts
finalNodeCount = len(G.nodes())
print 'There are ' + str(finalNodeCount) + ' nodes after removing nodes without published dates'
print 'There are ' + str(len(G.edges())) + ' edges after removing nodes without published dates'
print

# Print the oldest, newest, and range of dates
if (newestDate.month >= oldestDate.month):
    totalMonths = 12*(newestDate.year - oldestDate.year) + (newestDate.month - oldestDate.month)
else:
    totalMonths = 12*(newestDate.year - oldestDate.year) - (oldestDate.month - newestDate.month)

print 'Oldest published date: ' + str(oldestDate)
print 'Newest published date: ' + str(newestDate)
print 'Date range: ' + str(totalMonths) + ' months'
print

# Centrality
finalCentralities = nx.centrality.in_degree_centrality(G)
finalCentralNode = max(finalCentralities, key=finalCentralities.get)
print 'The most central node in the final graph is: ' + str(finalCentralNode) + ' with centrality ' + str(finalCentralities[finalCentralNode])
print 'Published on ' + G.node[finalCentralNode]['published']
print

# In degree: max, average, median
finalInDegrees = G.in_degree()
finalInDegreeMax = max(finalInDegrees, key=finalInDegrees.get)
print 'The highest in degree node in the final graph is: ' + str(finalInDegreeMax) + ' with in degree ' + str(finalInDegrees[finalInDegreeMax])
print 'Published on ' + G.node[finalInDegreeMax]['published']

averageFinalInDegree = sum(finalInDegrees.values())/float(len(finalInDegrees))
print 'The average in degree in the final graph is: ' + str(averageFinalInDegree)
print

sortedFinalInDegrees = sorted(finalInDegrees.values())
if not finalNodeCount % 2:
    medianFinalInDegree = (sortedFinalInDegrees[finalNodeCount / 2] + sortedFinalInDegrees[finalNodeCount / 2 - 1]) / 2.0
else: 
    medianFinalInDegree = sortedFinalInDegrees[finalNodeCount / 2]    
print 'The median in degree in the final graph is: ' + str(medianFinalInDegree)
print

# Find all nodes with at least 1 citation
allNodeInDegreeRates = {}
allInDegreeNodes = filter(lambda x: G.in_degree(x) > 0, G.nodes())
print str(len(allInDegreeNodes)) + ' papers are cited at least once'

# Graph the cumulative total of nodes with in degree <= x
allInDegreeValues = G.in_degree(allInDegreeNodes).values()
minInDegree = min(allInDegreeValues)
maxInDegree = max(allInDegreeValues)
cumulativeInDegreeNodeTotals = []

while (minInDegree <= maxInDegree):
    cumulativeInDegreeNodeTotals.append(len(filter(lambda x: x <= minInDegree, allInDegreeValues)))
    minInDegree = minInDegree + 1

cumulativePlot = subplot(111)
cumulativePlot.set_title("Cumulative total of nodes with in degree less than or equal to x")
cumulativePlot.set_xlabel("Log of In degree")
cumulativePlot.set_ylabel("Cumulative total of nodes with in degree <= x")
cumulativePlot.scatter(map(math.log, range(min(allInDegreeValues), maxInDegree + 1)), cumulativeInDegreeNodeTotals)
savefig("cumulativeInDegreeTotal.png")
close()

#######################################
#
# Step through the graph a month at a time
# Find cohort data
# Find nodes that were most central at some point
#
#######################################

firstOfMonth = date(newestDate.year, newestDate.month, 1) # Date to track which nodes to remove
G_temp = deepcopy(G)

mostCentralTimestepNodes = [] # List of the most central nodes at each timestep
cohortSize = [] # Size of each cohort
averageCohortOutDegree = [] # Average bibliography length of each cohort
while (len(G_temp.nodes()) > 0):
    
    # Grab the most central nodes by in degree
    centralities = nx.centrality.in_degree_centrality(G_temp)
    sortedCentralities = sorted(centralities, key=lambda x: centralities[x], reverse=True)
    i = 10
    if (len(sortedCentralities) >= i):
        while (i > 0):
            mostCentralTimestepNodes.append(sortedCentralities[i])
            i = i - 1
    else:
        for n in sortedCentralities:
            mostCentralTimestepNodes.append(n)

    # Find the nodes to remove
    removeNodes = filter(lambda x: makeDate(G_temp.node[x]['published']) > firstOfMonth, G_temp)

    # Find the size of this cohort
    cohortSize.append(len(removeNodes))
    
    # Find the average out degree of this cohort
    averageCohortOutDegree.append(sum(G.out_degree(removeNodes).values()) / float(len(removeNodes)))

    # Update the temporary graph and cutoff month
    G_temp.remove_nodes_from(removeNodes)
    if (firstOfMonth.month == 1):
        firstOfMonth = date(firstOfMonth.year - 1, 12, 1)
    else:
        firstOfMonth = date(firstOfMonth.year, firstOfMonth.month - 1, 1)

# Sort the cohort lists in chronological order
averageCohortOutDegree.reverse()
cohortSize.reverse();

# Plot the average cohort size
cohortSizePlot = subplot(111)
cohortSizePlot.set_title("Cohort size")
cohortSizePlot.set_xlabel("Month from the start")
cohortSizePlot.set_ylabel("Size of the cohort")
cohortSizePlot.scatter(range(1, len(cohortSize) + 1), cohortSize)
savefig("cohortSize.png")
close()

# Plot the average out degree by cohort
outDegreePlot = subplot(111)
outDegreePlot.set_title("Average out degree by month cohort")
outDegreePlot.set_xlabel("Month from the start")
outDegreePlot.set_ylabel("Average cohort out degree")
outDegreePlot.scatter(range(1, len(averageCohortOutDegree) + 1), averageCohortOutDegree)
savefig("cohortOutDegreeAverage.png")
close()

# Grab the unique list of most central nodes
mostCentralNodes = list(set(mostCentralTimestepNodes))

# Remove most central nodes with in degree < 15
# These were likely only most central because of small graph size
centralNodesToRemove = []
for m in mostCentralNodes:
    if (G.in_degree(m) < 15):
        centralNodesToRemove.append(m)
        
for r in centralNodesToRemove:
    mostCentralNodes.remove(r)

print 'We have ' + str(len(mostCentralNodes)) + ' different nodes that are most central at some point'
print 

#######################################
#
# Step through the graph a month at a time again
# Find the increase in in degree per timestep for each most central node
#
#######################################

centralNodeInDegreeRates = {}
firstOfMonth = date(newestDate.year, newestDate.month, 1)
G_temp = deepcopy(G)

while (len(G_temp.nodes()) > 0):
    
    # Find the nodes to remove
    removeNodes = filter(lambda x: makeDate(G_temp.node[x]['published']) > firstOfMonth, G_temp)
            
    # Find the in degree rates for the most central nodes
    for n in mostCentralNodes:
        
        # Initialize lists if necessary
        if (not n in centralNodeInDegreeRates.keys()):
            centralNodeInDegreeRates[n] = []
            
        # Add the in degree if the node exists
        if (n in G_temp.nodes()):
            centralNodeInDegreeRates[n].append(G_temp.in_degree(n))
        else:
            centralNodeInDegreeRates[n].append(0)

    # Update the temporary graph
    G_temp.remove_nodes_from(removeNodes)
    if (firstOfMonth.month == 1):
        firstOfMonth = date(firstOfMonth.year - 1, 12, 1)
    else:
        firstOfMonth = date(firstOfMonth.year, firstOfMonth.month - 1, 1)
    
# Get the rates of change in in degree for the most central nodes
# We only got the in degree at each step when iterating through the time steps
# Subtract previous total to get the change
for n in centralNodeInDegreeRates.keys():
    
    centralNodeInDegreeRates[n].reverse()
    previousCentralInDegree = 0
    tempDegreeRates = []
    
    for i in range(len(centralNodeInDegreeRates[n])):
        tempDegreeRates.append(centralNodeInDegreeRates[n][i] - previousCentralInDegree)        
        previousCentralInDegree = centralNodeInDegreeRates[n][i]
        
    centralNodeInDegreeRates[n] = deepcopy(tempDegreeRates)  

#################################
#
# Find the months to first citation
#
#################################

# Find months until first citation for each node with > 0 in links, average
monthsToFirstCitation = []
totalInDegreeFirstCitation = []
for n in allInDegreeNodes:

    #get all the neighbors
    allNeighbors = G.predecessors(n)
    firstCitationDate = date.max
    for neighbor in allNeighbors:
        neighborDate = makeDate(G.node[neighbor]['published'])
        if (neighborDate < firstCitationDate):
            firstCitationDate = neighborDate

    publishedStr = G.node[n]['published'].split('-')
    publishedDate = date(int(publishedStr[0]), int(publishedStr[1]), 1)
    
    if (firstCitationDate.month < publishedDate.month):
        x = (12 * (firstCitationDate.year - publishedDate.year) - (publishedDate.month - firstCitationDate.month))
    else:
        x = (12 * (firstCitationDate.year - publishedDate.year) + (firstCitationDate.month - publishedDate.month))

    # Consider only nodes referenced by nodes that are newer
    # Some nodes appear to be cited by papers that were published before them
    if (x > 0):
        monthsToFirstCitation.append(x)
        totalInDegreeFirstCitation.append(G.in_degree(n))    

# Print out some properties of months to first citation for all nodes with > 0 citations
averageFirstMonthToCitation = sum(monthsToFirstCitation) / float(len(monthsToFirstCitation))
print 'Average first month to citation for all papers: ' + str(averageFirstMonthToCitation)
print 'Range all: ' + str(max(monthsToFirstCitation) - min(monthsToFirstCitation))
print 'Max all: ' + str(max(monthsToFirstCitation))
print 'Min all: ' + str(min(monthsToFirstCitation))
print 'std deviation all: ' + str(np.std(monthsToFirstCitation))
print

firstCitationTotalCitations = subplot(111)
firstCitationTotalCitations.set_title("Total citations by months to first citation")
firstCitationTotalCitations.set_xlabel("Months to first citation")
firstCitationTotalCitations.set_ylabel("Total citations")
firstCitationTotalCitations.scatter(monthsToFirstCitation, totalInDegreeFirstCitation)
savefig("firstCitationByTotalCitations_all.png")
close()

firstCitationTotalCitations = subplot(111)
firstCitationTotalCitations.set_title("Total citations by months to first citation")
firstCitationTotalCitations.set_xlabel("Months to first citation")
firstCitationTotalCitations.set_ylabel("Log of Total citations")
firstCitationTotalCitations.scatter(monthsToFirstCitation, map(math.log, totalInDegreeFirstCitation))
savefig("firstCitationByTotalCitations_all_log.png")
close()

# Find months until first citation for each most central, average
# Copy centralNodeInDegreeRates to shift to published date
monthsToFirstCitationCentral = {}
adjustedMonthsToFirstCitationCentral = deepcopy(centralNodeInDegreeRates)

for n in mostCentralNodes:
    indegrees = centralNodeInDegreeRates[n]
    monthCounter = 0
    for i in indegrees:
        if (i > 0):
            break
        monthCounter = monthCounter + 1
        
    years = monthCounter / 12
    months = monthCounter - 12*years
    if (oldestDate.month + months > 12):
        firstCitationDateCentral = date(oldestDate.year + years + 1, (oldestDate.month + months) % 12, 1)
    else:
        firstCitationDateCentral = date(oldestDate.year + years, oldestDate.month + months, 1)
        
    publishedStr = G.node[n]['published'].split('-')
    publishedDate = date(int(publishedStr[0]), int(publishedStr[1]), 1)
    
    if (firstCitationDateCentral.month < publishedDate.month):
        monthsToFirstCitationCentral[n] = 12 * (firstCitationDateCentral.year - publishedDate.year) - (publishedDate.month - firstCitationDateCentral.month)
    else:
        monthsToFirstCitationCentral[n] = 12 * (firstCitationDateCentral.year - publishedDate.year) + (firstCitationDateCentral.month - publishedDate.month)

# Print out some properties of months to first citation for most central nodes
averageFirstMonthToCitationCentral = sum(monthsToFirstCitationCentral.values()) / float(len(monthsToFirstCitationCentral))
print 'Average: ' + str(averageFirstMonthToCitationCentral)
print 'Range: ' + str(max(monthsToFirstCitationCentral.values()) - min(monthsToFirstCitationCentral.values()))
print 'std deviation: ' + str(np.std(monthsToFirstCitationCentral.values()))
print 'Range: ' + str(max(monthsToFirstCitationCentral.values()) - min(monthsToFirstCitationCentral.values()))
print 'Max: ' + str(max(monthsToFirstCitationCentral.values()))
print 'Min: ' + str(min(monthsToFirstCitationCentral.values()))
print

# Shift months to first citation to publication date
adjustedCentralNodeInDegreeRates = deepcopy(centralNodeInDegreeRates)
for n in adjustedCentralNodeInDegreeRates:

    # Find how many months to first citation
    tempMonths = deepcopy(adjustedCentralNodeInDegreeRates[n])
    pubDateStr = G.node[n]['published'].split('-')
    pubDate = date(int(pubDateStr[0]), int(pubDateStr[1]), 1)

    if (pubDate.month < oldestDate.month):
        monthDiff = 12 * (pubDate.year - oldestDate.year) - (oldestDate.month - pubDate.month)
    else:
        monthDiff = 12 * (pubDate.year - oldestDate.year) + (pubDate.month - oldestDate.month)

    # Shift rates to start with published date
    # Pad out the back side of these lists so we can graph them
    d = deque(tempMonths)
    while (monthDiff > 0):
        t = d.popleft()
        d.append(t)
        monthDiff = monthDiff - 1

    tempMonths = list(d)
    adjustedCentralNodeInDegreeRates[n] = deepcopy(tempMonths)

# Find the average in degree from published date for most central nodes
averageAdjustedCentralNodeInDegreeRates = []
for i in range(0, len(adjustedCentralNodeInDegreeRates[adjustedCentralNodeInDegreeRates.keys()[0]])):
    ave = 0
    for n in adjustedCentralNodeInDegreeRates:
        ave = ave + adjustedCentralNodeInDegreeRates[n][i]
    averageAdjustedCentralNodeInDegreeRates.append(ave / float(len(mostCentralNodes)))

    inDegreePlot = subplot(111)
    inDegreePlot.set_title("Average In Degree acceleration curve for \nnodes that were most central at some point")
    inDegreePlot.set_xlabel("Month from publication date")
    inDegreePlot.set_ylabel("Average new in degrees")
    scatter(range(1, len(averageAdjustedCentralNodeInDegreeRates) + 1), averageAdjustedCentralNodeInDegreeRates)
    savefig("AverageInDegreeCurve.png")
    close()
    
    

    
