
secNames = [("Agriculture, Hunting, Forestry and Fishing",1),
			("Mining and Quarrying",2),
			("Food, Beverages and Tobacco",3),
			("Textiles and Leather",4),
			("Wood and Products of Wood and Cork",5),
			("Pulp, Paper, Paper , Printing and Publishing",6),
			("Coke, Refined Petroleum and Nuclear Fuel",7),
			("Chemicals, Rubber, and Plastics",8),
			("Other Non-Metallic Mineral",9),
			("Basic Metals and Fabricated Metal",10),
			("Machinery, Nec",11),
			("Electrical and Optical Equipment",12),
			("Transport Equipment",13),
			("Manufacturing, Nec; Recycling",14),
			("Electricity, Gas and Water Supply",15),
			("Construction",16),
			("Sales, Hotels, and Restaurants",17),
			("Inland Transport",18),
			("Water Transport",19),
			("Air Transport",20),
			("Post and Telecommunications",21),
			("Financial Intermediation",22),
			("Real Estate",23),
			("Other Public Expenditure",24),
			("Other Private Expenditure",25)]

secNames = DataFrame(sec=map(x->x[2],secNames),name=map(x->x[1],secNames))
save(dirs.data*"secNames.csv",secNames)
