import mysql.connector

def insert_mutation(name, chromosome, version, position, deleted_seq, inserted_seq, clinical_significance, user_inserted):
    try:
        conn = mysql.connector.connect(
            host='172.20.128.1',
            user='root',
            password='password',
            database='test'
        )
        cursor = conn.cursor()

        # Insert query      v - table name
        query = """
        INSERT INTO mutation_sorted ( 
           NAME, Chromosome, Version, Position,
            deleted_seq, inserted_seq, CLINICAL_SIGNIFICANCE, User_inserted
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
        """

        values = (
            name, chromosome, version, position,
            deleted_seq, inserted_seq, clinical_significance, user_inserted
        )
        
        cursor.execute(query, values)
        conn.commit()
        print(f"Mutation inserted successfully.")
    
    except mysql.connector.Error as err:
        print(f"Error: {err}")
    
    finally:
        if 'conn' in locals() and conn.is_connected():
            cursor.close()
            conn.close()


insert_mutation('BRCA1', 1, 0, 1, 'A', 'T', 'Pathogenic', 0)