USE rtx_demo_db;
DROP TABLE IF EXISTS Msmts; 
DROP TABLE IF EXISTS Channels; 
# DROP TABLE IF EXISTS Distros;

CREATE TABLE Channels (  
			id		INTEGER AUTO_INCREMENT, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE,
			r_lim	DOUBLE,
			info	VARCHAR(255),
			
			PRIMARY KEY (id)
			 
			); 
CREATE TABLE Msmts (
			id		INTEGER AUTO_INCREMENT, 
			time	DATETIME NOT NULL, 
			value	DOUBLE NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);

SELECT * FROM Channels