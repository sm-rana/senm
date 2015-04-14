set search_path to net1_scada,public;
DROP TABLE IF EXISTS msmts; 
DROP TABLE IF EXISTS channels; 
-- DROP TABLE IF EXISTS Distros;

CREATE TABLE channels (  
			id		SERIAL, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE PRECISION,
			r_lim	DOUBLE PRECISION,
			stde	DOUBLE PRECISION,
			info	VARCHAR(255),
			
			PRIMARY KEY (id)
			 
			); 
create index on channels (id);

CREATE TABLE msmts (
			id		SERIAL, 
			time	TIMESTAMP NOT NULL, 
			value	DOUBLE PRECISION NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);
create index on msmts (id);

-- Net1 "Scada"  system
-- Control components
INSERT INTO channels (net_id_str, msmt_t) VALUES 
	('PU9', 'C') ;

-- Online sensors - tank levels and limits
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim) VALUES
	('R9', 'L', 0, 0),  -- reservior
	('T2', 'L', 950, 1000) ;


-- Online sensors - flowrates
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim, stde) VALUES
	('PU9', 'F', -1000, 1000, 100), 

	('P21', 'Q', -2000, 2000, 20),
	 ('P31', 'Q', -2000, 2000, 15);

-- Online sensors - pressures
INSERT INTO channels (net_id_str, msmt_t, stde) VALUES
    ('10', 'B', 10),

	('12', 'P', 10), 
    ('23', 'P', 10),
	('32', 'P', 10);

-- online sensor - real-time water meters

INSERT INTO channels (net_id_str, msmt_t) VALUES
    ('31', 'D');



	
