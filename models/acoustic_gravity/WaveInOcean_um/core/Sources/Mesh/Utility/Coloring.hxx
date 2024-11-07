
namespace OndoMathX {
   
//! Creates the graph of connections by the vertices of the mesh
inline void getConnectivityMatrix(const Mesh & mesh, vector<vector<Index> > & connectivityMatrix)
{
    
    //First we make a list that contains for each Coord the number of elements that touch this Coord
    std::map<Index, std::vector<Index> > verticesElements;
    
 
    
    for (Index i=0;i<mesh.getNumElements();++i)
    {
        const Index numV = mesh.getElement(i).getNumVertices();
        
        for (Index j=0;j<numV;++j)
        {
            Index index_v  = mesh.getElement(i).getVertex(j).getIndex();
            
            verticesElements[index_v].push_back(i);
        }
    }
    
    //Second step : we (over) fill the connections matrix
    connectivityMatrix.resize(mesh.getNumElements());
    
    for (Index i=0;i<mesh.getNumElements();++i)
    {
        const Index numV = mesh.getElement(i).getNumVertices();
        
        for (Index j=0;j<numV;++j)
        {
            Index index_v  = mesh.getElement(i).getVertex(j).getIndex();
            const Index numVV = verticesElements[index_v].size();
            
            for (Index k=0;k<numVV;++k)
            {
                Index elem = verticesElements[index_v][k];
                if (elem != i) connectivityMatrix[i].push_back(elem);
            }
            
        }
    }
    
    //Final step remove the elements written more than one time
    for (Index i=0;i<mesh.getNumElements();++i)
    {
        sort(connectivityMatrix[i].begin(), connectivityMatrix[i].end());
        auto it = unique (connectivityMatrix[i].begin(), connectivityMatrix[i].end());
        connectivityMatrix[i].resize(std::distance(connectivityMatrix[i].begin(),it));
    }
}


//Get the coloring of the elements (using the connectivity of the highest dimension element)
inline void getColoring(const Mesh & mesh, std::vector<Index> & coloring, Index & numColors)
{
    Index dimension = mesh.getDimElement();
    Index no_color = mesh.getNumElements()+1;
    
    //construct the connectivity matrix
    std::vector<std::vector<Index> > connectivityMatrix;
    getConnectivityMatrix(mesh,connectivityMatrix);
    
    //Resize the output array filled with a sufficiently large Index indicating that no color has been used
    coloring.resize(connectivityMatrix.size(),no_color);
    
    //List of color description
    //first : number elements for this color
    //second : number of the color
    std::vector<std::pair<Index,Index> > size_color;
    
    //Initialization of the array with the first color (color number = 0)
    size_color.push_back(std::make_pair(0,0));
    
    //Number of color used
    numColors = 0;
    
    //We recover the first element of the proper dimension if the list of elements
    Index i_start = 0;
    while(mesh.getElement(i_start).getDim() != dimension) i_start++;
    
    //We compute the index of the element of the adequate dimension with the highest degree of connectivity ( degree of connectivity = number of neighbourg ).
    for (Index i = 0;i< connectivityMatrix.size();++i)
    {
        if (mesh.getElement(i).getDim() == dimension  && (connectivityMatrix[i].size()>connectivityMatrix[i_start].size()))
        {
            i_start = i;
        }
    }


 

    //Construct a queue initialize with the starting element
    std::queue<Index> F;
    F.push(i_start);
    
    //Main loop
    while(!F.empty())
    {
        //Selecting the element to deal with
        Index i_now = F.front();
        F.pop();
        
        Index c_tmp = 0;
        Index color = 0;
        bool test = true;
        
        //Check if the element has a color
        if (coloring[i_now] == no_color)
        {
            //We determine a color (we possibly create a new color)
            while(test && (c_tmp <= numColors))
            {
                test=false;
                color = size_color[c_tmp].second;
                
                //We find a suitable color to put the element in
                for (Index k=0;(k<connectivityMatrix[i_now].size()) && (!test);k++)
                {
                    //We test if the neighbour has the same color
                    if (coloring[connectivityMatrix[i_now][k]] == color)
                    {
                        c_tmp++;
                        test=true;
                    }
                }
            }
            
            //If there is no color available : we add a new one
            if (c_tmp > numColors)
            {
                numColors++;
                
                color = numColors;
                c_tmp = numColors;
                
                size_color.push_back(std::make_pair(0,numColors));
            }
            
            size_color[c_tmp].first++;
            sort(size_color.begin(),size_color.end());
            
            coloring[i_now]=color;
            
            //We add the neigbours in the queue
            for (Index k=0;k<connectivityMatrix[i_now].size();k++)
            {
                if ((coloring[connectivityMatrix[i_now][k]] == no_color)
                    && (mesh.getElement(connectivityMatrix[i_now][k]).getDim() == dimension))
                    F.push(connectivityMatrix[i_now][k]);
            }
        }
    }

    //Increment the number of colors for consistency
    numColors++;
}
     
    
} // namespace OndoMathX




